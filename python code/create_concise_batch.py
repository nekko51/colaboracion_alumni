import os
import shutil
import re
 
def natural_sort_key(text):
    #split into text & numbers
    parts = re.split(r'(\d+)', text)
    result = []
    for part in parts:
        if part.isdigit():
            #Convert to integer
            result.append(int(part))
        else:
            #Lowercase
            result.append(part.lower())
    return result
 
def find_latest_batch(metropolis_path):
    #find last batch ran (last one in alphabetical order)
    try:
        #list all the contents of directory
        contents = os.listdir(metropolis_path)
        batch_dirs = []
        #only keep if directory
        for item in contents:
            full_path = os.path.join(metropolis_path, item)
            if os.path.isdir(full_path):
                batch_dirs.append(item)
 
        if not batch_dirs:
            #if none, error message
            print("No batch folders were found.")
            return None
        
        #sort by name and keep last one
        batch_dirs.sort()
        latest_batch = batch_dirs[-1]
        
        return os.path.join(metropolis_path, latest_batch)
    except FileNotFoundError:
        print(f"Error: Could not find the folder '{metropolis_path}'")
        return None
 
def process_run_file(file_lines):
    #reads run_X.txt and extracts the required info (seed, seed energy, best seq)
 
    original_seed = None
    original_energy = None
    
    best_sequence = None
    best_energy = float('inf')
    best_hamming_dist = None
 
    #parse original seed
    for line in file_lines:
        if "seed:" in line and original_seed is None:
            match = re.search(r"seed: (\S+)", line)
            if match:
                original_seed = match.group(1)
        
        if "total energy:" in line and original_energy is None:
            #first total energy is that of the seed
            match = re.search(r"total energy:\s+([\d.-]+)", line)
            if match:
                original_energy = match.group(1)
        
        if original_seed and original_energy:
            break
 
    #search for lowest energy sequence of the entire run
    current_energy = None
    current_hamming_dist = None
    for line in file_lines:
        if line.strip().startswith("delta_e:"):
            #total energy and hamming dist
            energy_match = re.search(r"total energy:\s+([\d.-]+)", line)
            if energy_match:
                try:
                    current_energy = float(energy_match.group(1))
                except:
                    current_energy = None
            
            hamming_match = re.search(r"hamming distance to original:\s+([\d.-]+)", line)
            if hamming_match:
                current_hamming_dist = hamming_match.group(1)
            else:
                current_hamming_dist = None
        
        elif line.strip().startswith("resulting chain:"):
            if current_energy is not None:
                sequence_match = re.search(r"resulting chain: (\S+)", line)
                if sequence_match:
                    sequence = sequence_match.group(1)
                    
                    #look for lowest energy
                    if current_energy < best_energy:
                        best_energy = current_energy
                        best_sequence = sequence
                        best_hamming_dist = current_hamming_dist
                
                current_energy = None
                current_hamming_dist = None
 
    if best_energy == float('inf'):
        best_energy = None
 
    return original_seed, original_energy, best_sequence, best_energy, best_hamming_dist
 
 
def main():
    #script assumes it's in the 'python code' folder
    base_dir = os.path.join(os.path.dirname(__file__), '..')
    results_dir = os.path.join(base_dir, 'results')
    metropolis_dir = os.path.join(results_dir, 'metropolis')
 
 
    print("Searching for latest sim batch folder...")
    latest_batch = find_latest_batch(metropolis_dir)
 
    if not latest_batch:
        print("\nError: No batch found in 'results/metropolis'.")
        return
 
    batch_name = os.path.basename(latest_batch)
    print(f"Batch found: {batch_name}")
 
    #create summary folder
    output_dir = os.path.join(results_dir, 'metropolis_concise', f'concise_batch_{batch_name}')
    os.makedirs(output_dir, exist_ok=True)
 
    #copy betas file
    betas_source_path = os.path.join(latest_batch, 'betas.txt')
    if os.path.exists(betas_source_path):
        shutil.copy(betas_source_path, output_dir)
 
    #process results
    summary_path = os.path.join(output_dir, "concise_summary.txt")
    ultra_summary_path = os.path.join(output_dir, "ultra_concise_summary.txt")
    print(f"Creating summary file '{summary_path}'...")
    print(f"Creating ultra-concise summary file '{ultra_summary_path}'...")
    with open(summary_path, 'w') as summary_file, open(ultra_summary_path, 'w') as ultra_summary_file:
        
        #look for all seq folders
        seq_folders = []
        for item in os.listdir(latest_batch):
            if item.startswith('seq_') and os.path.isdir(os.path.join(latest_batch, item)):
                seq_folders.append(item)
        seq_folders.sort(key=natural_sort_key)
 
        for folder in seq_folders:
            seq_path = os.path.join(latest_batch, folder)
            summary_file.write(f"==================== {folder.upper()} ====================\n\n")
 
            run_files = []
            for f in os.listdir(seq_path):
                if f.startswith('run_') and f.endswith('.txt'):
                    run_files.append(f)
            run_files.sort(key=natural_sort_key)            
            
            if not run_files:
                summary_file.write("No 'run' files were found for this sequence.\n\n\n")
                continue
 
            #seed is same for all runs of a given seq
            with open(os.path.join(seq_path, run_files[0]), 'r') as f:
                lineas = f.readlines()
            original_seed, original_energy, _, _, _ = process_run_file(lineas)
 
            summary_file.write(f"Original seed:\n")
            summary_file.write(f"{original_seed}\n")
            summary_file.write(f"Original energy: {original_energy}\n\n")
 
            absolute_best_sequence = None
            absolute_best_energy = float('inf')
            absolute_best_hamming = None
            absolute_best_run_name = None

            #look for best result
            for run_name in run_files:
                with open(os.path.join(seq_path, run_name), 'r') as f:
                    lineas = f.readlines()
                
                #best step of THAT specific run
                _, _, best_sequence, best_energy, best_hamming = process_run_file(lineas)
 
                summary_file.write(f"Best step from '{run_name}':\n")
                if best_sequence:
                    summary_file.write(f"{best_sequence}\n")
                    summary_file.write(f"Energy: {best_energy}\n")
                    summary_file.write(f"Hamming distance to original: {best_hamming}\n\n")

                    if best_energy is not None  and  best_energy < absolute_best_energy:
                        absolute_best_energy = best_energy
                        absolute_best_sequence = best_sequence
                        absolute_best_hamming = best_hamming
                        absolute_best_run_name = run_name
                else:
                    summary_file.write("Could not find a valid result in this run.\n\n")
 
            if absolute_best_sequence:
                ultra_summary_file.write(f"Best result for {folder.upper()} (from '{absolute_best_run_name}'):\n\n")
                ultra_summary_file.write(f"Seed: {original_seed}\n")
                ultra_summary_file.write(f"Best: {absolute_best_sequence}\n")
                ultra_summary_file.write(f"Energy: {absolute_best_energy}\n")
                ultra_summary_file.write(f"Hamming distance to original: {absolute_best_hamming}\n\n\n")
            else:
                ultra_summary_file.write(f"Could not find any valid result for {folder.upper()}.\n\n\n")

            print(f"Processed folder: {folder}")
 
    print("\nDone!")

if __name__ == "__main__":
    main()