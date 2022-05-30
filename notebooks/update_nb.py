import os
import subprocess


for file in os.listdir(os.getcwd()):
    if file.endswith(".ipynb"):
        print(f'updating {file}')
        proc = subprocess.run(
            ['python', 'nb_convert.py', file],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL
        )