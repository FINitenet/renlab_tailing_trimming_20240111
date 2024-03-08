#!/usr/bin/python3
# -*- coding:utf-8 -*-
# @FileName  :tailing_trimming_profile.py
# @Time      :2024/03/08 15:00:47
# @Author    :Yuchen@rlab
# @Description: This script is used to profile the tailing trimming results of miRNA-seq data.

import os
import csv
import subprocess
import datetime
import argparse
import numpy as np
import pandas as pd
from Bio.Seq import reverse_complement

class Profile_maker:
    def __init__(self, meta_file, input_file, output_path, sample_name):
        self.metafile = meta_file
        self.script = os.path.join(os.path.dirname(__file__), "4_163-profile.pl")
        self.sample_name = sample_name
        self.output_path = output_path
        self.input_file = input_file # remapping result: merged-alignment-sorted.sam
        self.output_file = os.path.join(self.output_path, sample_name + ".txt")
        
    def dyc_163_profile(self):
        print(f"[{datetime.datetime.now()}] Start to generate the miRNA profile!")
        with open(self.output_file, "w") as f:
            profile_cmd = f"perl self.script {self.metafile} {self.input_file}"
            subprocess.run(profile_cmd, check=True, shell=True, stdout=f)
        print(f"[{datetime.datetime.now()}] Profile generated!")

    def split_163_file(self):
        print(f"[{datetime.datetime.now()}] Start to split the 163 profile file!")
        with open(self.output_file, "r") as f:
            lines = f.readlines()
        matrix = [line.strip().split('\t') for line in lines]
        groups = [matrix[i:i+12] for i in range(0, len(matrix), 12)]
        for group in groups:
            last_line = group.pop()
            newgroup = [row[3:] for row in group[3:]]
            newgroup.insert(0, last_line) 
            group.insert(0, last_line) 
            new_filename = newgroup[0][0].split()[0].replace("*", "_star") + ".txt"
            file = os.path.join(self.output_path, self.sample_name, new_filename)
            with open(file, "w") as f:
                f.write('\t\n'.join(['\t'.join(row) for row in group]) + '\t\n')

    def martix_summary(self):
        n_lines = len(open(self.output_file, 'r').readlines())
        n_groups = int(n_lines/12)
        res = []
        for i in range(n_groups):
            temp_res = []
            matrix_data = np.loadtxt(self.output_file, skiprows=i*12, max_rows=11)
            tag_data = np.loadtxt(self.output_file, skiprows=i*12+11, max_rows=1, dtype=str)[0]
            temp_res.append(tag_data)
            sum_all = np.sum(matrix_data)
            temp_res.append(sum_all)
            temp_res.append(str(sum_all - np.sum(matrix_data[:,10])))
            temp_res.append(str(sum_all - np.sum(matrix_data[10,:])))
            # for j in range(10):
            for j in range(9,-1,-1):
                temp_res.append(str(np.sum(matrix_data[j,:])))
            res.append(temp_res)

        res = np.array(res)
        np.savetxt(os.path.join(self.output_path, self.sample_name + ".summary.txt"), res, fmt="%s", delimiter = "\t",
                    header = "ID\tSUM\ttrimming\ttailling\ttail_1\ttail_2\ttail_3\ttail_4\ttail_5\ttail_6\ttail_7\ttail_8\ttail_9\ttail_10")

    def data_process(self):
        print(f"[{datetime.datetime.now()}] Start to process the profile data!")
        self.dyc_163_profile()
        self.split_163_file()
        self.martix_summary()
        print(f"[{datetime.datetime.now()}] Profile data processed!")


class Profile_summary:
    def __init__(self, meta_file, input_file, output_path, sample_name):
        self.metafile = meta_file
        self.sample_name = sample_name
        self.input_file = input_file # remapping result: merged-alignment-sorted.sam
        self.output_path = output_path

    def get_mirna_5GMC_reads(self):
        id = []
        site = []
        with open(self.metafile, mode='r') as f:
            data = csv.reader(f, delimiter='\t')
            for d in data:
                id.append(d[0])
                site.append(d[2])
        f.close()

        out = []
        with open(self.input_file, mode='r') as f:
            data = csv.reader(f, delimiter='\t')
            for d in data:
                for i, s in zip(id, site):
                    if i.startswith(d[2]) and d[3] == s:
                        d[2] = i
                        if d[2] == '16':
                            d[9] = reverse_complement(d[9])
                        new_line = '\t'.join([d[0], d[2], d[3], d[9]]) + '\n' # type: ignore
                        out.append(new_line)
                        break
        f.close()
        
        output_file = os.path.join(self.output_path, self.sample_name + ".5GMC")
        header = "RNAME\tID\tStart\tSequence\n"
        out.insert(0, header)
        with open(output_file, 'w') as f:
            f.writelines(out)
        return output_file

    def get_mirna_5GMC_sequence(self, output_file):
        df = pd.read_csv(output_file, sep='\t')
        groups = df.groupby(df.columns[1])
        for name, group in groups:
            filename = os.path.join(self.output_path, "doc_seqlogo", self.sample_name, name.replace("*", "_star") + ".txt") # type: ignore
            group.iloc[:, 3].apply(lambda x: x.replace('T', 'U')).to_csv(filename, sep='\t', index=False, header=False)

        with open(self.metafile, 'r') as f:
            for line in f:
                file_name = line.strip().split("\t")[0].replace("*","_star")
                file_path = os.path.join(self.output_path, file_name + ".txt")
                if not os.path.exists(file_path):
                    with open(file_path, 'w') as f:
                        f.write('N'*30)
    
    def data_process(self):
        print(f"[{datetime.datetime.now()}] Start to summary the profile data!")
        output_file = self.get_mirna_5GMC_reads()
        self.get_mirna_5GMC_sequence(output_file)
        print(f"[{datetime.datetime.now()}] Profile data processed!")

##### Argument parsing / help message /version

parser = argparse.ArgumentParser(prog=os.path.basename(__file__))
parser.add_argument('-i', '--input', help='inputdir containing merged-alignment-sorted.sam', required=True)
parser.add_argument('-o', '--output', help='output path', required=True)
parser.add_argument('-f', help='miRNA meta file, contains miRNA start end length', 
                    default = os.path.join(os.path.dirname(__file__), "source/miRNA_start_sequence_length_new.txt"),required=True)
parser.add_argument('-v', '--version', action='version', version='%(prog)s 20240307')

args = parser.parse_args()

##### Script control .. go through each file and process

if __name__ == "__main__":
    # RUN PATH
    INPUTDIR = args.input
    OUTDIR = args.output
    METAFILE = args.f
    PROFILE_PATH = os.path.join(OUTDIR, "4_163.results")
    SUMMARY_PATH = os.path.join(OUTDIR, "5_GMC_analysis")

    INPUT_FILES = os.path.join(INPUTDIR, "merged-alignment-sorted.sam") 
    SAMPLE_NAME = os.path.basename(INPUTDIR)
    
    # Create the output directory
    if not os.path.exists(PROFILE_PATH):
        os.makedirs(PROFILE_PATH)
    if not os.path.exists(SUMMARY_PATH):
        os.makedirs(SUMMARY_PATH)
    
    # Profile the tailing trimming results
    profile_maker = Profile_maker(METAFILE, INPUT_FILES, PROFILE_PATH, SAMPLE_NAME)
    profile_maker.data_process()

    # Summary the tailing trimming results
    profile_summary = Profile_summary(METAFILE, INPUT_FILES, SUMMARY_PATH, SAMPLE_NAME)
    profile_summary.data_process()

    print(f"[{datetime.datetime.now()}] the {SAMPLE_NAME} profile data has been processed!")