import requests
import os
import subprocess
import shutil


def alphafold_downloader(pdb_id, path):
    pdb_id = pdb_id.rstrip()
    url = "https://alphafold.ebi.ac.uk/files/AF-"+pdb_id+"-F1-model_v2.pdb"
    print(url)
    r = requests.get(url, allow_redirects=True).content.decode("utf-8")
    with open(f"{path}/temp/{pdb_id}.pdb", "w") as f:
        f.write(r)
    return pdb_id


def FASTA_iterator(fasta_filename):
    actual_protein = None
    sequence = ""
    with open(fasta_filename) as f:
        for line in f.readlines():
            line = line.rstrip()
            if line.startswith(">"):
                if actual_protein is None:
                    actual_protein = line
                    continue
                return tuple([actual_protein, sequence])
                actual_protein = line
                sequence = ""
                continue
            sequence += line
        return tuple([actual_protein, sequence])


# PSI-BLAST
def get_psiblast_prot(sequence, database, path):
    """This function returns the result of the psi-blast REST API"""
    shutil.copyfile(f"{path}/blastscript/psiblast.py",
                    f"{path}/temp/psiblast.py")
    os.chdir('./temp/')
    p = subprocess.Popen(['python', f"{path}/temp/psiblast.py", '--email', 'pyflex@protonmail.com',
                          '--sequence', sequence, '--database', database], stdout=subprocess.PIPE)
    print("Psi-blast is running...")
    i = ""
    for line in p.stdout:
        i = line
    p.wait()
    print(p.returncode)
    os.chdir('../')
    print(os.getcwd())
    return i


def decode_job_id(bytedata):
    """This function returns the job ID as a string, because the result of the previously subprocess is always a bytes object"""
    job_id = bytedata.decode('UTF-8')
    print("JOB ID \n", job_id, "\n")
    job_id = job_id.split(" ")
    job_id = job_id[3].split(".")
    return job_id[0]


def best_target_candidate(file_with_candidates):
    with open(file_with_candidates) as f:
        file = f.readlines()
    candidate1 = file[0].split(":")
    return candidate1[1]
