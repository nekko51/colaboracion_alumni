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
    path_parser = re.compile(r'betas(\d+)_iterations(\d+)')

    with open(output_file_path, 'w') as outfile:
        # Escribir la cabecera del CSV
        # Primero, leemos el primer archivo para saber cuántas lambdas hay
        with open(lambda_files_found[0], 'r') as f:
            num_lambdas = len(f.read().strip().split('\t')) - 1
        
        header = "n_betas,n_iterations,energy," + ",".join([f"lambda_{i}" for i in range(num_lambdas)])
        outfile.write(header + "\n")

        # Procesar cada archivo encontrado
        for lambda_path in lambda_files_found:
            match = path_parser.search(lambda_path)
            if not match:
                print(f"AVISO: No se pudieron extraer los parámetros de la ruta: {lambda_path}")
                continue

            n_betas = int(match.group(1))
            n_iterations = int(match.group(2))

            with open(lambda_path, 'r') as infile:
                line = infile.read().strip()
                if not line:
                    print(f"AVISO: El archivo está vacío: {lambda_path}")
                    continue
                
                # Los valores están separados por tabuladores. Los reemplazamos por comas para el formato CSV.
                values = line.replace('\t', ',').rstrip(',')
                outfile.write(f"{n_betas},{n_iterations},{values}\n")

    print(f"\n¡Éxito! Se han agrupado {len(lambda_files_found)} archivos en '{output_file_path}'.")

if __name__ == "__main__":
    aggregate_svm_lambdas()
    aggregate_svm_summaries()
