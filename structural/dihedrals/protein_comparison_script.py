# PyMOL script to compare two aligned protein structures
# Calculates differences in CA distances and backbone angles (phi, psi, omega)

from pymol import cmd
import numpy as np
import math

def calc_ca_distance_diff(struct1, struct2, output_file="ca_distances.txt"):
    """
    Calculate differences in CA-CA distances between two structures
    """
    print(f"Calculating CA distance differences between {struct1} and {struct2}...")
    
    # Get CA coordinates for both structures
    coords1 = []
    coords2 = []
    
    # Extract CA coordinates
    cmd.iterate_state(1, f"{struct1} and name CA", "coords1.append([x,y,z])", space=locals())
    cmd.iterate_state(1, f"{struct2} and name CA", "coords2.append([x,y,z])", space=locals())
    
    coords1 = np.array(coords1)
    coords2 = np.array(coords2)
    
    if len(coords1) != len(coords2):
        print(f"Warning: Different number of CA atoms ({len(coords1)} vs {len(coords2)})")
        min_len = min(len(coords1), len(coords2))
        coords1 = coords1[:min_len]
        coords2 = coords2[:min_len]
    
    # Calculate pairwise distances for each structure
    n_residues = len(coords1)
    distances1 = np.zeros((n_residues, n_residues))
    distances2 = np.zeros((n_residues, n_residues))
    
    for i in range(n_residues):
        for j in range(n_residues):
            distances1[i,j] = np.linalg.norm(coords1[i] - coords1[j])
            distances2[i,j] = np.linalg.norm(coords2[i] - coords2[j])
    
    # Calculate differences
    distance_diff = distances2 - distances1
    
    # Save results
    with open(output_file, 'w') as f:
        f.write("Residue_i\tResidue_j\tDistance_Struct1\tDistance_Struct2\tDifference\n")
        for i in range(n_residues):
            for j in range(i+1, n_residues):  # Only upper triangle
                f.write(f"{i+1}\t{j+1}\t{distances1[i,j]:.3f}\t{distances2[i,j]:.3f}\t{distance_diff[i,j]:.3f}\n")
    
    print(f"CA distance differences saved to {output_file}")
    return distance_diff

def calc_angle_diff(struct1, struct2, output_file="angle_differences.txt"):
    """
    Calculate differences in phi, psi, and omega angles between two structures
    """
    print(f"Calculating backbone angle differences between {struct1} and {struct2}...")
    
    # Get residue information for both structures
    residues1 = []
    residues2 = []
    
    # Extract residue info (chain, residue number, residue name)
    cmd.iterate(f"{struct1} and name CA", "residues1.append((chain, resi, resn))", space=locals())
    cmd.iterate(f"{struct2} and name CA", "residues2.append((chain, resi, resn))", space=locals())
    
    results = []
    
    for i, (chain, resi, resn) in enumerate(residues1):
        try:
            # Calculate angles for structure 1
            phi1 = cmd.get_dihedral(f"{struct1} and chain {chain} and resi {int(resi)-1} and name C",
                                   f"{struct1} and chain {chain} and resi {resi} and name N",
                                   f"{struct1} and chain {chain} and resi {resi} and name CA",
                                   f"{struct1} and chain {chain} and resi {resi} and name C")
            
            psi1 = cmd.get_dihedral(f"{struct1} and chain {chain} and resi {resi} and name N",
                                   f"{struct1} and chain {chain} and resi {resi} and name CA",
                                   f"{struct1} and chain {chain} and resi {resi} and name C",
                                   f"{struct1} and chain {chain} and resi {int(resi)+1} and name N")
            
            omega1 = cmd.get_dihedral(f"{struct1} and chain {chain} and resi {int(resi)-1} and name CA",
                                     f"{struct1} and chain {chain} and resi {int(resi)-1} and name C",
                                     f"{struct1} and chain {chain} and resi {resi} and name N",
                                     f"{struct1} and chain {chain} and resi {resi} and name CA")
            
            # Calculate angles for structure 2
            phi2 = cmd.get_dihedral(f"{struct2} and chain {chain} and resi {int(resi)-1} and name C",
                                   f"{struct2} and chain {chain} and resi {resi} and name N",
                                   f"{struct2} and chain {chain} and resi {resi} and name CA",
                                   f"{struct2} and chain {chain} and resi {resi} and name C")
            
            psi2 = cmd.get_dihedral(f"{struct2} and chain {chain} and resi {resi} and name N",
                                   f"{struct2} and chain {chain} and resi {resi} and name CA",
                                   f"{struct2} and chain {chain} and resi {resi} and name C",
                                   f"{struct2} and chain {chain} and resi {int(resi)+1} and name N")
            
            omega2 = cmd.get_dihedral(f"{struct2} and chain {chain} and resi {int(resi)-1} and name CA",
                                     f"{struct2} and chain {chain} and resi {int(resi)-1} and name C",
                                     f"{struct2} and chain {chain} and resi {resi} and name N",
                                     f"{struct2} and chain {chain} and resi {resi} and name CA")
            
        except:
            # If angle calculation fails, set to None
            phi1 = phi2 = psi1 = psi2 = omega1 = omega2 = None
        
        # Calculate differences (handle circular nature of angles)
        def angle_diff(a1, a2):
            if a1 is None or a2 is None:
                return None
            diff = a2 - a1
            while diff > 180:
                diff -= 360
            while diff < -180:
                diff += 360
            return diff
        
        phi_diff = angle_diff(phi1, phi2)
        psi_diff = angle_diff(psi1, psi2)
        omega_diff = angle_diff(omega1, omega2)
        
        results.append({
            'chain': chain,
            'resi': resi,
            'resn': resn,
            'phi1': phi1, 'phi2': phi2, 'phi_diff': phi_diff,
            'psi1': psi1, 'psi2': psi2, 'psi_diff': psi_diff,
            'omega1': omega1, 'omega2': omega2, 'omega_diff': omega_diff
        })
    
    # Save results
    with open(output_file, 'w') as f:
        f.write("Chain\tResidue\tResName\tPhi1\tPhi2\tPhi_Diff\tPsi1\tPsi2\tPsi_Diff\tOmega1\tOmega2\tOmega_Diff\n")
        for r in results:
            f.write(f"{r['chain']}\t{r['resi']}\t{r['resn']}\t")
            
            # Handle None values gracefully
            phi_str = f"{r['phi1']:.2f}\t{r['phi2']:.2f}\t{r['phi_diff']:.2f}\t" if r['phi_diff'] is not None else "N/A\tN/A\tN/A\t"
            psi_str = f"{r['psi1']:.2f}\t{r['psi2']:.2f}\t{r['psi_diff']:.2f}\t" if r['psi_diff'] is not None else "N/A\tN/A\tN/A\t"
            omega_str = f"{r['omega1']:.2f}\t{r['omega2']:.2f}\t{r['omega_diff']:.2f}\n" if r['omega_diff'] is not None else "N/A\tN/A\tN/A\n"
            
            f.write(phi_str + psi_str + omega_str)
    
    print(f"Angle differences saved to {output_file}")
    return results

def analyze_structural_differences(struct1, struct2, output_prefix="analysis"):
    """
    Complete analysis of structural differences between two aligned proteins
    """
    print(f"\n=== ANALYZING STRUCTURAL DIFFERENCES ===")
    print(f"Structure 1: {struct1}")
    print(f"Structure 2: {struct2}")
    print("=" * 50)
    
    # Calculate CA distance differences
    distance_diff = calc_ca_distance_diff(struct1, struct2, f"{output_prefix}_ca_distances.txt")
    
    # Calculate angle differences
    angle_results = calc_angle_diff(struct1, struct2, f"{output_prefix}_angles.txt")
    
    # Print summary statistics
    print("\n=== SUMMARY STATISTICS ===")
    print(f"Average CA distance difference: {np.mean(np.abs(distance_diff)):.3f} Å")
    print(f"Max CA distance difference: {np.max(np.abs(distance_diff)):.3f} Å")
    print(f"RMS CA distance difference: {np.sqrt(np.mean(distance_diff**2)):.3f} Å")
    
    # Calculate angle statistics
    phi_diffs = [r['phi_diff'] for r in angle_results if r['phi_diff'] is not None]
    psi_diffs = [r['psi_diff'] for r in angle_results if r['psi_diff'] is not None]
    omega_diffs = [r['omega_diff'] for r in angle_results if r['omega_diff'] is not None]
    
    if phi_diffs:
        print(f"Average |phi| angle difference: {np.mean(np.abs(phi_diffs)):.2f}°")
        print(f"Max |phi| angle difference: {np.max(np.abs(phi_diffs)):.2f}°")
    
    if psi_diffs:
        print(f"Average |psi| angle difference: {np.mean(np.abs(psi_diffs)):.2f}°")
        print(f"Max |psi| angle difference: {np.max(np.abs(psi_diffs)):.2f}°")
    
    if omega_diffs:
        print(f"Average |omega| angle difference: {np.mean(np.abs(omega_diffs)):.2f}°")
        print(f"Max |omega| angle difference: {np.max(np.abs(omega_diffs)):.2f}°")
    
    print(f"\nOutput files:")
    print(f"- {output_prefix}_ca_distances.txt")
    print(f"- {output_prefix}_angles.txt")
    print("=" * 50)

# Example usage - MODIFY THESE NAMES TO MATCH YOUR STRUCTURES
def run_comparison():
    """
    Main function to run the comparison
    CHANGE THE STRUCTURE NAMES BELOW TO MATCH YOUR PYMOL OBJECTS
    """
    # Replace these with your actual PyMOL object names
    structure1_name = "4dkl"  # <-- CHANGE THIS
    structure2_name = "5c1m"  # <-- CHANGE THIS
    output_prefix = "my_analysis"   # <-- CHANGE THIS IF DESIRED
    
    # Run the analysis
    analyze_structural_differences(structure1_name, structure2_name, output_prefix)

# Uncomment the line below to run automatically when script is loaded
run_comparison()

print("Script loaded successfully!")
print("To run the analysis, either:")
print("1. Modify the structure names in run_comparison() and uncomment the last line, then reload")
print("2. Or call: analyze_structural_differences('your_struct1_name', 'your_struct2_name')")