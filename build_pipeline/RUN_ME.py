import subprocess
import os

scripts = ['protein_list_builder.py', 'link_to_uniprot.py', 'clean_data.py', 'link_with_fasta.py']

repo_dir = os.path.expanduser('~/HPP_TargetList/build_pipeline')

for script in scripts:
    script_path = os.path.join(repo_dir, script)
    print(f"Running {script_path}")
    try:
        subprocess.run(['python3', script_path], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error while running {script_path}: {e}")
