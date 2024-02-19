# Give a amino acid sequence, the script will output elemental composition, mono-isotopic & avg-mass
import numpy as np
import pandas as pd

# General script description...
print("This script will output the elemental composition, monoisotopic mass, & average mass for an amino acid sequence.")
print("\nAtomic weights & isotopic composition were pulled from the NIST database:")
print("https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses")

# seq_request = 'peptides'

aa_mass = pd.read_csv("aminoacidtable.csv", encoding="windows-1252")
element_isotopes = pd.read_csv("elementisotopes.csv", encoding="windows-1252")

# Function for elemental composition calculation
def seq_composition(sequence):
    '''function calculates monoisotopic mass for a sequence'''
    sequence = sequence.upper()
    aa_freq = {}
    # Create frequency dictionary of ammino acids
    for i in sequence:
        if i in aa_freq:
            aa_freq[i] += 1
        else:
            aa_freq[i] = 1

    element_sum = {'C':0,'H':0,'N':0,'O':0, 'S':0}
    merged = {}
    for key, value in aa_freq.items():

        count = int(value)
    
        aa_dict = aa_mass.query(f"Code == '{key}'")["Residue Dictionary"].to_dict()
        
        for key in aa_dict:
            
            updated_aa_dict = eval(aa_dict[key])
            
            for key in updated_aa_dict:
                updated_aa_dict[key] *= count
            
        for key in element_sum:
            if key in updated_aa_dict:
                element_sum[key] = element_sum[key] + updated_aa_dict[key]
            else:
                pass

    # Add remaining atoms for N- & C-terminus            
    element_sum["H"] += 2
    element_sum["O"] += 1

    return element_sum

# Function for monoisotopic mass calculation
def mono_mass(sequence):
    '''returns monoisotopic mass for amino acid sequence'''
    composition = seq_composition(sequence)
    
    # Iterate through elemental composition and multiply number by isotopic mass
    total_isomass = 0.0
    for key in composition:
        atom = key
        
        isotopic_mass = element_isotopes.query(f"Code == '{key}'")["Isotopic Mass"].to_dict()
        iso_mass = [float(x) for x in list(isotopic_mass.values())]

        for key in composition:
            if key == atom:
                mass = composition[f'{atom}'] * iso_mass[0]
                
                total_isomass = total_isomass + mass
            
    return "%.4f" % round(total_isomass,4)

# Function for average mass calculation
def average_mass(sequence):
    '''returns the average mass of a amino acid sequence'''
    composition = seq_composition(sequence)

    total_avgmass = 0.0

    for key, value in composition.items():

        atom_count = value


        isotope_dict = element_isotopes.query(f"Code == '{key}'")["Mass Abundance"].to_dict()

        avg_atom_mass = 0
        for key, value in isotope_dict.items():

            extracted_iso = eval(value)


            for key, value in extracted_iso.items():
                avg_atom_mass += (key * value * atom_count)
                

            total_avgmass += avg_atom_mass

    return "%.4f" % round(total_avgmass, 4)


# Request peptide or protein sequence
seq_request = input("Enter amino acid sequence: ")

# Limited amino acid lettes (for version 1)
request_lim = ['A','R','N','D','C','E','Q','G','H','O','I','L','K','M','F','P','U','S','T','W','Y','V']

# Test if sequence contains letters outside of list
for i in seq_request.upper():
    if i in request_lim:
        request_code = True
    else:
        request_code = False
    
if request_code == True:
# Print returned information for requested amino acid sequence
    print("\nThe elemental composition for the sequence is:", seq_composition(seq_request))
    print("\nThe monoisotopic mass (Da):", mono_mass(seq_request))
    print("\nThe average mass (Da):", average_mass(seq_request))
# Print error message
elif request_code == False:
    print("\nSequence contains a non-coding letter.")


input("\nPress the Enter key to continue... ") 
