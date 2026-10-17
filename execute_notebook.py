"""Execute the report's read/analysis cells without rerunning the MD trajectories."""
import os
from pathlib import Path

root = Path(__file__).resolve().parent
os.chdir(root)
for key, directory in {
    'JUPYTER_RUNTIME_DIR': root / '.tmp/jupyter-runtime',
    'IPYTHONDIR': root / '.tmp/ipython',
    'JUPYTER_CONFIG_DIR': root / '.tmp/jupyter-config',
    'JUPYTER_DATA_DIR': root / '.venv/share/jupyter',
    'TEMP': root / '.tmp',
    'TMP': root / '.tmp',
}.items():
    directory.mkdir(parents=True, exist_ok=True)
    os.environ[key] = str(directory)
os.environ['JUPYTER_PATH'] = str(root / '.venv/share/jupyter')

import nbformat
from nbclient import NotebookClient

notebook_path = root / '铜热膨胀研究.ipynb'
notebook = nbformat.read(notebook_path, as_version=4)
client = NotebookClient(notebook, timeout=120, kernel_name='copper-md',
                        resources={'metadata': {'path': str(root)}})
client.execute()
errors = [output for cell in notebook.cells if cell.cell_type == 'code'
          for output in cell.get('outputs', []) if output.output_type == 'error']
if errors:
    raise RuntimeError(errors)
nbformat.validate(notebook)
temporary = root / '.tmp/report_executed.ipynb'
nbformat.write(notebook, temporary)
temporary.replace(notebook_path)
print('Notebook executed successfully, no cell errors')
