import subprocess

from glob import glob

command = "fastp -i {0} -o {1} -I {2} -O {3} --detect_adapter_for_pe"

files = glob('./data/*R1*')
for file in files:
    base = file.split('_R1')[0]
    r1 = [file, f'{base}_trimmed_R1.fastq.gz']
    r2 = [file.replace('R1', 'R2'), f'{base}_trimmed_R2.fastq.gz']
    com = command.format(*r1, *r2)
    subprocess.run(com, shell=True)


