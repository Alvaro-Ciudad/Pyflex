#!/usr/bin/env python

import time
import os
import sys
import shutil
from interface_ui import *
from qt_material import apply_stylesheet
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from pyflex import downloader, protein_treatment, GNMBuild, outputfun
import warnings
from Bio import BiopythonWarning


class MainWindow(QtWidgets.QMainWindow, Ui_MainWindow):

    def __init__(self, *args, **kwargs):
        QtWidgets.QMainWindow.__init__(self, *args, **kwargs)
        self.setupUi(self)
        self.pushButton.clicked.connect(self.btnstart_clicked)
        self.progressBar.hide()
        self.label_img.hide()

    def get_fasta_sequence(self):
        fasta_sequence = self.plainTextEdit_1.toPlainText()
        return fasta_sequence

    def get_fasta_file_path(self):
        fasta_file = self.plainTextEdit_2.toPlainText()
        return fasta_file
    
    def get_output_name(self):
        out_name = self.plainTextEdit_3.toPlainText()
        return out_name        

    def btnstart_clicked(self):
        fasta_sequence = self.get_fasta_sequence()
        fasta_file_path = self.get_fasta_file_path()
        output_name = self.get_output_name()
        self.worker = WorkerThread(fasta_sequence, fasta_file_path, output_name)
        self.worker.start()
        self.worker.worker_complete.connect(self.worker_finished)
        self.worker.update_progress.connect(self.evt_update_progress)
    
    def worker_finished(self, emp):
        QMessageBox.information(self, "Done!", "PyFlex completed\n\n{}".format(emp["message"]))
        self.label_pyflex.hide()
        self.pushButton.hide()
        self.progressBar.hide()
        self.label_submit.hide()
        self.pixmap = QPixmap(f"/outputs/{self.get_output_name()}_graph.png")
        self.label_img.setPixmap(self.pixmap)
        self.label_img.show()

    def evt_update_progress(self, val):
        self.progressBar.show()
        self.label_submit.setText("Your job is succesfully submitted!")
        self.label_input.hide()
        self.label_output.hide()
        self.label_in1.hide()
        self.label_in2.hide()
        self.label_out1.hide()
        self.plainTextEdit_1.hide()
        self.plainTextEdit_2.hide()
        self.plainTextEdit_3.hide()
        self.pushButton.hide()
        self.progressBar.setValue(val)


class WorkerThread(QThread):
    update_progress = pyqtSignal(int)
    worker_complete = pyqtSignal(dict)

    def __init__(self, fasta_seq, fasta_path, out_name):
        QThread.__init__(self)
        self.path = os.path.dirname(os.path.abspath(__file__))
        self.fasta_seq = fasta_seq
        self.fasta_path = fasta_path
        self.out_name = out_name


    def run(self):
        x = 0
        while x < 100:
            self.update_progress.emit(x)
    
            if self.fasta_path and self.fasta_seq:
                raise ValueError("Fill only one of the two input options.\n")

            else:
                if os.path.isfile(self.fasta_path):
                    id, seq = downloader.FASTA_iterator(self.fasta_path)
                    x += 10
                    self.update_progress.emit(x)

                    job_id, ret_code = downloader.get_psiblast_prot(
                        seq, "uniprotkb_swissprot", self.path)
                    if ret_code != 0:
                        shutil.rmtree('temp')
                        raise ValueError("Problems with Psi-Blast, restart the program.")
                    job = downloader.decode_job_id(job_id)
                    x += 10
                    self.update_progress.emit(x)

                    best_candidate = protein_treatment.best_target_candidate(
                        f"{self.path}/temp/{job}.preselected_ids.txt")
                    x += 10
                    self.update_progress.emit(x)

                    id = downloader.alphafold_downloader(best_candidate, self.path)
                    x += 10
                    self.update_progress.emit(x)

                    scores = GNMBuild.obtain_MSFs(f"{self.path}/temp/{id}.pdb")
                    x += 10
                    self.update_progress.emit(x)

                    job_id_pdb, ret_code = downloader.get_psiblast_prot(seq, "pdb", self.path)
                    if ret_code != 0:
                        shutil.rmtree('temp')
                        raise ValueError("Problems with Psi-Blast, restart the program.")
                    job_pdb = downloader.decode_job_id(job_id_pdb)
                    x += 10
                    self.update_progress.emit(x)

                    # Protein treatment
                    files = protein_treatment.list_of_pdbfiles_preprocessor(
                        f"{self.path}/temp/{job_pdb}.ids.txt")
                    os.chdir('./temp/')
                    x += 10
                    self.update_progress.emit(x)

                    protein_treatment.pdb_download_files(files)
                    homo = protein_treatment.split_files_by_chain(files)
                    homo_object = protein_treatment.Homologues(homo)
                    x += 10
                    self.update_progress.emit(x)

                    # Output 
                    scaled_scores = GNMBuild.B_factor_scaling(homo_object, scores)
                    outputfun.write_pdb_CA(
                        best_candidate.rstrip(), scores, "simulated", self.path, self.out_name)
                    outputfun.write_pdb_CA(
                        best_candidate.rstrip(), scaled_scores, "scaled", self.path, self.out_name)
                    x += 10
                    self.update_progress.emit(x)
                    outputfun.plot_b_factors(scores, self.path, self.out_name)
                    os.chdir('../')
                    x += 10
                    self.update_progress.emit(x)
                    sys.stdout.write(
                        "Succesfully finished the job, output is in /outputs.\n")

                if len(self.fasta_seq) > 10:
                    job_id, ret_code = downloader.get_psiblast_prot(
                        self.fasta_seq, "uniprotkb_swissprot", self.path)
                    if ret_code != 0:
                        shutil.rmtree('temp')
                        raise ValueError("Problems with Psi-Blast, restart the program.")
                    
                    job = downloader.decode_job_id(job_id)
                    x += 10
                    self.update_progress.emit(x)

                    best_candidate = protein_treatment.best_target_candidate(
                    f"{self.path}/temp/{job}.preselected_ids.txt")
                    x += 10
                    self.update_progress.emit(x)

                    id = downloader.alphafold_downloader(best_candidate, self.path)
                    x += 10
                    self.update_progress.emit(x)

                    scores = GNMBuild.obtain_MSFs(f"{self.path}/temp/{id}.pdb")
                    x += 10
                    self.update_progress.emit(x)

                    job_id_pdb, ret_code = downloader.get_psiblast_prot(
                        self.fasta_seq, "pdb", self.path)
                    if ret_code != 0:
                        shutil.rmtree('temp')
                        raise ValueError("Problems with Psi-Blast, restart the program.")
                    job_pdb = downloader.decode_job_id(job_id_pdb)
                    x += 10
                    self.update_progress.emit(x)

                    files = protein_treatment.list_of_pdbfiles_preprocessor(
                        f"{self.path}/temp/{job_pdb}.ids.txt")
                    os.chdir('./temp/')
                    x += 10
                    self.update_progress.emit(x)

                    protein_treatment.pdb_download_files(files)
                    homo = protein_treatment.split_files_by_chain(files)
                    homo_object = protein_treatment.Homologues(homo)
                    x += 10
                    self.update_progress.emit(x)

                    scaled_scores = GNMBuild.B_factor_scaling(homo_object, scores)
                    x += 10
                    self.update_progress.emit(x)

                    outputfun.write_pdb_CA(
                        best_candidate.rstrip(), scores, "simulated", self.path, self.out_name)
                    x += 10
                    self.update_progress.emit(x)
                    outputfun.write_pdb_CA(
                        best_candidate.rstrip(), scaled_scores, "scaled", self.path, self.out_name)

                    outputfun.plot_b_factors(scores, self.path, self.out_name)
                    os.chdir('../')
                    x += 10
                    self.update_progress.emit(x)

                    sys.stdout.write(
                        "Succesfully finished the job, output is in /outputs.\n")
                    

                else:
                    raise ValueError("Incorrect input, please introduce a sequence in FASTA format.\n")
        
        self.worker_complete.emit({"message":"Succesfully finished the job, output is in outputs directory."})

        


if __name__ == "__main__":

    try:
        os.mkdir("temp")
    except FileExistsError:
        shutil.rmtree('temp')
        os.mkdir("temp")

    warnings.simplefilter('ignore', BiopythonWarning)

    app = QtWidgets.QApplication([])
    apply_stylesheet(app, theme='dark_teal.xml')
    window = MainWindow()
    window.show()
    app.exec_()

    shutil.rmtree('temp')









