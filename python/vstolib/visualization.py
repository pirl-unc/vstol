import subprocess
from importlib import resources


def visualize(
        rscript_path: str,
        tsv_file: str,
        output_pdf_file: str
):
    """
    Visualize a TSV file.

    Args:
        rscript_path        :   Rscript path.
        tsv_file            :   VSTOL TSV file.
        output_pdf_file     :   Output PDF file.
    """
    with resources.path('vstolib.resources.scripts', 'visualize.R') as script_path:
        command = [rscript_path,
                   str(script_path),
                   '--tsv-file', tsv_file,
                   '--output-pdf-file', output_pdf_file]
        try:
            result = subprocess.run(command, check=True, capture_output=True, text=True)
            print("R script output:", result.stdout)
        except subprocess.CalledProcessError as e:
            print(f"Error running R script: {e.stderr}")
