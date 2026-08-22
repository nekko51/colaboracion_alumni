import os
import re


def aggregate_svm_summaries():
    """
    Busca todos los archivos 'summary.txt' dentro del directorio 'results/svm/',
    y los agrupa en un único archivo de texto para facilitar el análisis.
    """
    # La ruta base asume que el script está en 'python code/'
    base_dir = os.path.join(os.path.dirname(__file__), '..')
    svm_results_dir = os.path.join(base_dir, 'results', 'svm', 'annealing')
    output_file_path = os.path.join(svm_results_dir, 'aggregated_svm_summary.txt')

    print(f"Buscando resúmenes en: {svm_results_dir}")
    print(f"El resultado se guardará en: {output_file_path}")

    summary_files_found = []
    # os.walk() recorre el árbol de directorios de arriba hacia abajo
    for dirpath, _, filenames in os.walk(svm_results_dir):
        if 'summary.txt' in filenames:
            summary_files_found.append(os.path.join(dirpath, 'summary.txt'))

    if not summary_files_found:
        print("No se encontraron archivos 'summary.txt'.")
        return

    # Ordena los archivos para un resultado consistente
    summary_files_found.sort()

    with open(output_file_path, 'w') as outfile:
        for summary_path in summary_files_found:
            # Obtiene una ruta relativa para que el encabezado sea más limpio
            relative_path = os.path.relpath(summary_path, base_dir)
            print(f"Agregando: {relative_path}")

            # Escribe un encabezado para identificar el origen de los datos
            outfile.write("=" * 80 + "\n")
            outfile.write(f"Contenido de: {relative_path}\n")
            outfile.write("=" * 80 + "\n\n")

            # Lee el contenido del summary.txt y lo escribe en el archivo de salida
            with open(summary_path, 'r') as infile:
                outfile.write(infile.read())
                outfile.write("\n\n")

    print(f"\n¡Éxito! Se han agrupado {len(summary_files_found)} archivos en '{output_file_path}'.")

def aggregate_svm_lambdas():
    """
    Busca todos los archivos 'result-lambdas.txt' dentro del directorio 'results/svm/annealing',
    extrae los parámetros de la ruta y agrupa los datos en un único archivo CSV.
    """
    # La ruta base asume que el script está en 'python code/'
    base_dir = os.path.join(os.path.dirname(__file__), '..')
    annealing_results_dir = os.path.join(base_dir, 'results', 'svm', 'annealing')
    output_file_path = os.path.join(annealing_results_dir, 'aggregated_svm_lambdas.csv')

    print(f"Buscando resultados lambda en: {annealing_results_dir}")
    print(f"El resultado se guardará en: {output_file_path}")

    lambda_files_found = []
    # os.walk() recorre el árbol de directorios
    for dirpath, _, filenames in os.walk(annealing_results_dir):
        if 'result-lambdas.txt' in filenames:
            lambda_files_found.append(os.path.join(dirpath, 'result-lambdas.txt'))

    if not lambda_files_found:
        print("No se encontraron archivos 'result-lambdas.txt'.")
        return

    # Ordena los archivos para un resultado consistente
    lambda_files_found.sort()

    # Expresión regular para extraer n_betas e iteraciones de la ruta
    # Ejemplo: .../betas00100_iterations01000/...
    path_parser_c_k = re.compile(r'betas(\d+)_iterations(\d+)_c-limit([0-9.eE+-]+)_kernel(\w)')
    path_parser_c = re.compile(r'betas(\d+)_iterations(\d+)_c-limit([0-9.eE+-]+)') # without kernel
    path_parser = re.compile(r'betas(\d+)_iterations(\d+)') # Fallback for older format

    with open(output_file_path, 'w') as outfile:
        # Escribir la cabecera del CSV
        # Primero, leemos el primer archivo para saber cuántas lambdas hay
        with open(lambda_files_found[0], 'r') as f:
            num_lambdas = len(f.read().strip().split('\t')) - 1
        
        header = "kernel,c_limit,n_betas,n_iterations,energy," + ",".join([f"lambda_{i}" for i in range(num_lambdas)])
        outfile.write(header + "\n")

        # Procesar cada archivo encontrado
        for lambda_path in lambda_files_found:
            kernel = 'g' # default
            match = path_parser_c_k.search(lambda_path)
            if match:
                n_betas = int(match.group(1))
                n_iterations = int(match.group(2))
                c_limit = float(match.group(3))
                kernel = match.group(4)
            else:
                match = path_parser_c.search(lambda_path)
                if match:
                    n_betas = int(match.group(1))
                    n_iterations = int(match.group(2))
                    c_limit = float(match.group(3))
                else:
                    match = path_parser.search(lambda_path)
                    if match:
                        n_betas = int(match.group(1))
                        n_iterations = int(match.group(2))
                        c_limit = 1.0 # Default C-limit for old format
                    else:
                        print(f"AVISO: No se pudieron extraer los parámetros de la ruta: {lambda_path}")
                        continue

            with open(lambda_path, 'r') as infile:
                line = infile.read().strip()
                if not line:
                    print(f"AVISO: El archivo está vacío: {lambda_path}")
                    continue
                
                # Los valores están separados por tabuladores. Los reemplazamos por comas para el formato CSV.
                values = line.replace('\t', ',').rstrip(',')
                outfile.write(f"{kernel},{c_limit},{n_betas},{n_iterations},{values}\n")

    print(f"\n¡Éxito! Se han agrupado {len(lambda_files_found)} archivos en '{output_file_path}'.")

def aggregate_decision_function_tests():
    """
    Busca todos los archivos 'test-results.txt' en 'results/svm/decision-function-tests',
    y agrupa los resultados en un único archivo de texto tabulado.
    """
    # La ruta base asume que el script está en 'python code/'
    base_dir = os.path.join(os.path.dirname(__file__), '..')
    tests_dir = os.path.join(base_dir, 'results', 'svm', 'decision-function-tests')    
    output_validation_path = os.path.join(tests_dir, 'aggregated_validation_results.txt')
    output_test_path = os.path.join(tests_dir, 'aggregated_test_results.txt')

    print(f"Buscando resultados de la función de decisión en: {tests_dir}")
    print(f"Los resultados de validación (667 filas) se guardarán en: {output_validation_path}")
    print(f"Los resultados de test (668 filas) se guardarán en: {output_test_path}")

    # Constantes para el número de filas esperado
    VALIDATION_ROWS = 667
    TEST_ROWS = 668

    test_files_found = []
    # os.walk() para encontrar todos los 'test-results.txt'
    for dirpath, _, filenames in os.walk(tests_dir):
        if 'test-results.txt' in filenames:
            test_files_found.append(os.path.join(dirpath, 'test-results.txt'))

    if not test_files_found:
        print("No se encontraron archivos 'test-results.txt'.")
        return

    test_files_found.sort()

    # Estructuras de datos para validación y test
    validation_tags, validation_results = [], {}
    test_tags, test_results = [], {}

    # Contadores para saber si ya hemos guardado los tags
    is_first_validation = True
    is_first_test = True

    # Contador para los nombres de columna simplificados, p.ej. {'c': 1, 'g': 2}
    kernel_counters = {}

    # Procesar cada archivo encontrado
    for i, test_path in enumerate(test_files_found):
        # Generar una cabecera de columna única y simplificada como 'c-001'.
        relative_dir = os.path.dirname(os.path.relpath(test_path, tests_dir))
        
        try:
            # Extrae el directorio base, p.ej. 'kernel-c' o 'kernel-w_C-0.90_W-0.20'
            kernel_dir = relative_dir.split(os.sep)[0]
            match = re.search(r'kernel-(\w)(?:_C-([0-9.]+)_W-([0-9.]+))?', kernel_dir)
            
            if match:
                kernel_char = match.group(1)
                if kernel_char == 'w' and match.group(2) and match.group(3):
                    center = match.group(2)
                    width = match.group(3)
                    header_base = f"w_C{center}_W{width}"
                else:
                    header_base = kernel_char
            else:
                header_base = kernel_dir.split('-')[-1] # Fallback

            # Incrementa el contador para esta configuración y formatea la cabecera
            kernel_counters[header_base] = kernel_counters.get(header_base, 0) + 1
            column_header = f"{header_base}-{kernel_counters[header_base]:03d}"
            
        except Exception as e:
            print(f"AVISO: No se pudo extraer el kernel de la ruta: {relative_dir}. Usando cabecera genérica.")
            column_header = f"unknown-{i+1:03d}"
        current_values = []
        try:
            # Primero, leemos todas las líneas para contar las filas
            with open(test_path, 'r') as infile:
                lines = infile.readlines()
            
            num_rows = len(lines)
            
            # Decidir si es un archivo de validación, de test, o ninguno
            if num_rows == VALIDATION_ROWS:
                target_tags, target_results, is_first = validation_tags, validation_results, is_first_validation
            elif num_rows == TEST_ROWS:
                target_tags, target_results, is_first = test_tags, test_results, is_first_test
            else:
                print(f"AVISO: El archivo '{test_path}' tiene {num_rows} filas, no coincide con validación ({VALIDATION_ROWS}) ni con test ({TEST_ROWS}). Saltando.")
                continue

            # Procesar las líneas del archivo
            for line in lines:
                line_parts = line.strip().split('\t')
                
                if not line_parts or not line_parts[0]: continue

                if len(line_parts) >= 3: # Formato nuevo con 3 columnas: [tag, dec_func, tag*dec_func]
                    tag, value = line_parts[0], line_parts[2] # Usamos la tercera columna pre-calculada
                    if is_first: target_tags.append(tag)
                    current_values.append(value)
                elif len(line_parts) == 2: # Formato con 2 columnas: [tag, dec_func]
                    tag, dec_func_val = line_parts[0], line_parts[1]
                    # Calculamos el valor nosotros mismos
                    value = float(dec_func_val) * int(tag)
                    if is_first: target_tags.append(tag)
                    # Lo guardamos como string para mantener la consistencia
                    current_values.append(f"{value:.15e}")
                else: # Formato antiguo: [value]
                    value = line_parts[0]
                    current_values.append(value)
            
            target_results[column_header] = current_values

            # Actualizar los flags de "primer archivo"
            if num_rows == VALIDATION_ROWS: is_first_validation = False
            if num_rows == TEST_ROWS: is_first_test = False

        except Exception as e:
            print(f"Error procesando el archivo {test_path}: {e}")

    # Escribir los archivos de salida
    def write_aggregated_file(output_path, tags, results_data):
        if not results_data:
            print(f"No se encontraron datos para agrupar en '{output_path}'.")
            return
        
        with open(output_path, 'w') as outfile:
            headers = ["tag"] + list(results_data.keys())
            outfile.write("\t".join(headers) + "\n")

            for i in range(len(tags)):
                row_data = [tags[i]]
                for h in results_data.keys():
                    row_data.append(results_data[h][i])
                outfile.write("\t".join(row_data) + "\n")
        print(f"¡Éxito! Se han agrupado {len(results_data)} tests en '{output_path}'.")

    write_aggregated_file(output_validation_path, validation_tags, validation_results)
    write_aggregated_file(output_test_path, test_tags, test_results)


def calculate_quadratic_certainty(filepath):
    """
    Calculates the normalized Quadratic Certainty Score (Brier-based accuracy)
    from a 'test-results.txt' file, measuring both correctness and margin confidence.
    
    Metric formulation:
      B = 1 - (1/n) * sum_{i=1}^n (1 - sigma(s_i * f(x_i)))^2
    
    Returns a dictionary with accuracy, F1 scores, and quadratic certainty in [0, 1].
    """
    import math
    total = 0
    sq_err_sum = 0.0
    correct = 0
    tp = fp = fn = tn = 0

    with open(filepath, 'r') as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            tag = int(parts[0])
            if len(parts) >= 3:
                cert = float(parts[2]) # s * f
            elif len(parts) == 2:
                cert = tag * float(parts[1])
            else:
                continue

            # Sigmoid probability of the ground-truth class
            prob_true = 1.0 / (1.0 + math.exp(-min(max(cert, -500), 500)))
            sq_err_sum += (1.0 - prob_true) ** 2
            total += 1

            if cert > 0:
                correct += 1
                if tag == 1: tp += 1
                else: tn += 1
            else:
                if tag == 1: fn += 1
                else: fp += 1

    brier_acc = 1.0 - (sq_err_sum / total) if total > 0 else 0.0
    acc = correct / total if total > 0 else 0.0
    f1_h = 2 * tp / (2 * tp + fp + fn) if (2 * tp + fp + fn) > 0 else 0.0
    f1_m = 2 * tn / (2 * tn + fp + fn) if (2 * tn + fp + fn) > 0 else 0.0

    return {
        'filepath': filepath,
        'total': total,
        'accuracy': acc,
        'f1_human': f1_h,
        'f1_murine': f1_m,
        'quadratic_certainty': brier_acc,
        'tp': tp, 'tn': tn, 'fp': fp, 'fn': fn
    }


if __name__ == "__main__":
    aggregate_svm_lambdas()
    aggregate_svm_summaries()
    aggregate_decision_function_tests()
