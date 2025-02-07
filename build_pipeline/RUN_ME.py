import subprocess
import os

scripts = ['protein_list_builder.py', 'link_to_uniprot.py', 'manual_fix_list_maker.py','clean_entries.py', 'update_PE.py', 'peptideAtlas_link.py', 'link_with_fasta.py']

repo_dir = os.path.expanduser('~/HPP_TargetList/build_pipeline')

for script in scripts:
    script_path = os.path.join(repo_dir, script)
    print(f"Running {script_path}")
    try:
        subprocess.run(['python3', script_path], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error while running {script_path}: {e}")

