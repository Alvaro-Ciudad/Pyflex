import os
import argparse
from pyflex import downloader, protein_treatment, GNMBuild, outputfun
import warnings
from Bio import BiopythonWarning
import shutil


def parserfunc():
    parser = argparse.ArgumentParser(
        description="This program does blablablabla")

    requiredNamed = parser.add_argument_group("Required arguments")

    requiredNamed.add_argument('-i', '--input',
                               dest="infile",
                               action="store",
                               default=None,
                               required=True,
                               help="Input FASTA file/Uniprot Id")

    parser.add_argument('-o', '--output',
                        dest="outfile",
                        action="store",
                        default="output",
                        help="Output file")

    parser.add_argument('-g', '--graph',
                        dest="outgraph",
                        action="store",
                        default=None,
                        help="Graphical result")

    args = parser.parse_args()
    return args


def main():
    os.mkdir("temp")
    warnings.simplefilter('ignore', BiopythonWarning)
    arg = parserfunc()
    path = os.path.dirname(os.path.abspath(__file__))
    if os.path.isfile(arg.infile):
        id, seq = downloader.FASTA_iterator(arg.infile)
        job_id = downloader.get_psiblast_prot(seq, "uniprotkb_swissprot", path)
        print(job_id)
        job = downloader.decode_job_id(job_id)
        #job = "psiblast-R20220404-143845-0932-33667906-p2m"
        print(job)
        best_candidate = protein_treatment.best_target_candidate(
            f"{path}/temp/{job}.preselected_ids.txt")
        print("best", best_candidate)

        id = downloader.alphafold_downloader(best_candidate, path)
        print(id)
        scores = GNMBuild.obtain_MSFs(f"{path}/temp/{id}.pdb")
        job_id_pdb = downloader.get_psiblast_prot(seq, "pdb", path)
        job_pdb = downloader.decode_job_id(job_id_pdb)
        #job_pdb = "psiblast-R20220404-150732-0697-3296736-p2m"
        files = protein_treatment.list_of_pdbfiles_preprocessor(
            f"{path}/temp/{job_pdb}.ids.txt")
        print(files)
        os.chdir('./temp/')
        _ = protein_treatment.pdb_download_files(files)
        print(_)
        homo = protein_treatment.split_files_by_chain(files)
        print(homo)
        homo_object = protein_treatment.Homologues(homo)
        scaled_scores = GNMBuild.B_factor_scaling(homo_object, scores)
        print(scaled_scores)
        _ = outputfun.write_pdb_CA(
            best_candidate.rstrip(), scores, "simulated", path, arg.outfile)
        _ = outputfun.write_pdb_CA(
            best_candidate.rstrip(), scaled_scores, "scaled", path, arg.outfile)
        _ = outputfun.plot_b_factors(scores, path, arg.outfile)
        os.chdir('../')
        print(os.getcwd())
    elif len(arg.infile) > 10:
        job_id = downloader.get_psiblast_prot(
            arg.infile, "uniprotkb_swissprot", path)
        print(job_id)

        job = downloader.decode_job_id(job_id)
        print(job)

        best_candidate = protein_treatment.best_target_candidate(
            f"{path}/temp/{job}.preselected_ids.txt")

        id = downloader.alphafold_downloader(best_candidate, path)
        print(id)

        scores = GNMBuild.obtain_MSFs(f"{path}/temp/{id}.pdb")
        print(scores)
        job_id_pdb = downloader.get_psiblast_prot(arg.infile, "pdb", path)

    else:
        id = downloader.alphafold_downloader(arg.infile, path)
        print(id)

        scores = GNMBuild.obtain_MSFs(f"{path}/temp/{id}.pdb")

        job_id_pdb = downloader.get_psiblast_prot(arg.infile, "pdb", path)
    shutil.rmtree('temp')


if __name__ == '__main__':
    main()
