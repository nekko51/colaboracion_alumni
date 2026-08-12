from pathlib import Path

def build_dataset():
    # Obtiene el directorio donde está este script (python code)
    script_dir = Path(__file__).resolve().parent
    # Sube un nivel al directorio raíz del proyecto
    base_dir = script_dir.parent
    
    path_human = base_dir / 'seqs' / '00human.txt'
    path_mouse = base_dir / 'seqs' / '00mouse.txt'
    path_output = base_dir / 'seqs' / 'allchains.txt'

    # Crea el directorio si no existiera
    path_output.parent.mkdir(parents=True, exist_ok=True)

    with open(path_output, 'w') as f_out:
        if path_human.exists():
            with open(path_human, 'r') as f_human:
                for line in f_human:
                    seq = line.strip()
                    if seq:
                        f_out.write(f"1\t{seq}\n")

        if path_mouse.exists():
            with open(path_mouse, 'r') as f_mouse:
                for line in f_mouse:
                    seq = line.strip()
                    if seq:
                        f_out.write(f"0\t{seq}\n")

import random
from pathlib import Path

def count_labels(split_name, lines):
    human = sum(1 for line in lines if line.startswith('1'))
    mouse = sum(1 for line in lines if line.startswith('0'))
    print(f"{split_name} -> Humanas: {human}, Murinas: {mouse}")

def split_dataset(learn_perc, validation_perc):
    script_dir = Path(__file__).resolve().parent
    seqs_dir = script_dir.parent / 'seqs'
    
    path_all = seqs_dir / 'allchains.txt'
    path_learn = seqs_dir / 'learn.txt'
    path_val = seqs_dir / 'validation.txt'
    path_test = seqs_dir / 'test.txt'

    with open(path_all, 'r') as f:
        lines = f.readlines()

    lines = [line for line in lines if line.strip()]
    
    random.seed() 
    random.shuffle(lines)
    
    total_lines = len(lines)
    learn_idx = int(total_lines * (learn_perc / 100.0))
    val_idx = learn_idx + int(total_lines * (validation_perc / 100.0))
    
    learn_lines = lines[:learn_idx]
    val_lines = lines[learn_idx:val_idx]
    test_lines = lines[val_idx:]
    
    with open(path_learn, 'w') as f:
        f.writelines(learn_lines)
        
    with open(path_val, 'w') as f:
        f.writelines(val_lines)
        
    with open(path_test, 'w') as f:
        f.writelines(test_lines)

    count_labels("Learn", learn_lines)
    count_labels("Validation", val_lines)
    count_labels("Test", test_lines)

if __name__ == '__main__':
    split_dataset(70, 15)