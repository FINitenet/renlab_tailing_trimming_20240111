#!/usr/bin/python3
# -*- coding:utf-8 -*-
# @FileName  :tailing_trmming_preprocess.py
# @Time      :2024/03/07 15:18:47
# @Author    :Yuchen@rlab
# @Description: This script is used to preprocess the raw data of sRNA-seq, including trimming the adapter and low quality reads, and remapping the reads to the reference genome and trsnoRNA.

import os
import glob
import pysam
import gzip
import subprocess
import datetime
import argparse
import logging
from Bio import SeqIO
from Bio.Seq import Seq

class PreProcess:
    def __init__(self, input_file, output_path, adapter, sample_name, flag):
        self.input_file = input_file
        self.output_file = os.path.join(output_path, sample_name + "_trimmed.fasta")
        self.output_path = output_path 
        self.adapter = adapter
        self.sample_name = sample_name
        self.flag = flag
        self.log_file = os.path.join(output_path, "log.txt")

    def MeanPhredScore(self, trimmed_quality):
        score_sum = sum(ord(q) - 33 for q in trimmed_quality)
        mean_score = round(score_sum / len(trimmed_quality))
        return mean_score

    def extracted_mirna_with_quality_score(self):
        records = []
        phred_scores = [0] * 100
        with pysam.FastxFile(self.input_file, threads=24) as f:
            for record in f:
                sequence = record.sequence
                quality = record.quality
                if sequence.count(self.adapter) > 0:
                    trimmed_sequence = sequence.split(self.adapter, 1)[0]
                    trimmed_quality = quality[:len(trimmed_sequence)]
                else:
                    trimmed_sequence = sequence
                    trimmed_quality = quality
                
                if 30 >= len(trimmed_sequence) >= 12:
                    mean_score = self.MeanPhredScore(trimmed_quality)
                    phred_scores[int(mean_score)] += 1
                    if mean_score >= 20:
                        record_id = f"{self.sample_name}-{len(records)}-score:{mean_score}"
                        records.append((record_id, trimmed_sequence))
        for idx, count in enumerate(phred_scores):
            if count > 0:
                print(f"Score {idx}: {count} sequences")

        with open(self.output_file, "w") as output:
            SeqIO.write((SeqIO.SeqRecord(Seq(sequence), id=record_id, description="") for record_id, sequence in records), output, "fasta")
    
    def modify_fasta_header(self):
        output_records = []
        with gzip.open(self.input_file, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                record.id = f"{self.sample_name}-{len(output_records)}-score:37"
                record.description = ""
                output_records.append(record)
        with open(self.output_file, "w") as f:
            SeqIO.write(output_records, f, "fasta")
    
    def modify_fastq_header(self):
        output_records = []
        with gzip.open(self.output_path + '/' + self.sample_name + "_trimmed.fq.gz", "rt") as f:
            for record in SeqIO.parse(f, "fastq"):
                record.id = f"{self.sample_name}-{len(output_records)}-score:37"
                record.description = ""
                output_records.append(record)
        with open(self.output_file, "w") as f:
            SeqIO.write(output_records, f, "fasta")


    def data_process(self):
        log_file = open(self.log_file, "a")
        print(f"[{datetime.datetime.now()}] Start to trim the adapter and low quality reads!")
        fq_extensions = [".fastq", ".fq", ".fastq.gz", ".fq.gz"]
        fa_extensions = [".fasta", ".fa", ".fasta.gz", ".fa.gz"]
        if self.input_file.endswith(tuple(fq_extensions)) and self.flag == False:
            self.extracted_mirna_with_quality_score()
        elif self.input_file.endswith(tuple(fa_extensions)) and self.flag == False:
            self.modify_fasta_header()
        elif self.flag == True:
            trim_cmd = "trim_galore --fastqc --fastqc_args '-t 16 --nogroup' --gzip -q 20 --length 12 --max_length 30 --trim-n --basename {} --no_report_file -a {}  -j 8 -o {} {}".format(self.sample_name, self.adapter, self.output_path, self.input_file)
            subprocess.run(trim_cmd, check=True, shell=True, stdin = log_file, stderr= log_file)
            self.modify_fastq_header()
        else:
            print("Please check your input file format!")
            exit(1)
        log_file.close()

    def get_output_file(self):
        return self.output_file
    
class ReMapping:
    def __init__(self, output_path, reference, trsnorna, sample_name):
        self.input_file = preprocess.get_output_file()
        self.output_path = os.path.join(output_path, sample_name)
        self.reference = reference
        self.trsnoRNA = trsnorna
        self.sample_name = sample_name
        self.unmapped2genome = os.path.join(self.output_path, "unmapped-to-genome.fasta")
        self.unmapped2trsnorna = os.path.join(self.output_path, "clean.fasta")
        self.log_file = os.path.join(self.output_path, "remapping.log")

    def reformatting_fasta(self):
        print(f"[{datetime.datetime.now()}] Reformatting fasta header")
        output_records = []
        with open(self.input_file, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                record.id = record.id.replace("-score:37", f"-{record.seq}--0")
                record.description = ""
                output_records.append(record)
        with open(self.unmapped2genome, "w") as f:
            SeqIO.write(output_records, f, "fasta")

    def remapping(self):
        log_file = open(self.log_file, "w")
        print(f"Mapping whole reads to trsnoRNA and unmapped reads are saved to {self.output_path}/clean.fasta")
        mapping_cmd = "bowtie -p 24 -v 0 -S -a -x {} -f {} --un {}/clean.fasta {}/map2trsnoRNA.sam".format(self.trsnoRNA, self.unmapped2genome, self.output_path, self.output_path)
        subprocess.run(mapping_cmd, check=True, shell=True, stderr=log_file, stdout=log_file)
        print(f"Mapping whole reads to hairpin20 and unmapped reads are saved to {self.output_path}/unmapped-0.fasta")
        mapping_cmd = "bowtie -p 24 -v 0 -S -a -x {} -f {} --un {}/unmapped-0.fasta {}/mapped-0.sam".format(self.reference, self.unmapped2trsnorna, self.output_path, self.output_path)
        subprocess.run(mapping_cmd, check=True, shell=True, stderr=log_file, stdout=log_file)

        for i in range(1, 11):
            k = i - 1
            output_records = []
            print(f"Trimming {i} bp for each unmapped reads")
            with open(f"{self.output_path}/unmapped-{k}.fasta", "r") as in_file:
                for record in SeqIO.parse(in_file, "fasta"):
                    record.seq = record.seq[:-1]
                    if len(record.seq) >= 12:
                        record.id = record.id.split("--")[0] + f"--{i}"
                        output_records.append(record)

            with open(f"{self.output_path}/trimmed-{i}.fasta", "w") as out_file:
                SeqIO.write(output_records, out_file, "fasta")
            
            print(f"Mapping trimmed {i} bp reads to hairpin20 and unmapped reads are saved to {self.output_path}/unmapped-{i}.fasta")
            remapping_cmd = "bowtie -p 24 -v 0 -S -a --norc --no-unal -x {} -f {}/trimmed-{}.fasta --un {}/unmapped-{}.fasta {}/mapped-{}.sam".format(self.reference, self.output_path, i, self.output_path, i, self.output_path, i)
            subprocess.run(remapping_cmd, check=True, shell=True, stderr=log_file, stdout=log_file)
        log_file.close()
        

    def merge_result(self):
        for i in range(0, 11):
            print(f"Processing: {self.output_path}/mapped-{i}.sam")
            with pysam.AlignmentFile(f"{self.output_path}/mapped-{i}.sam", "r", threads=24) as in_file:
                with pysam.AlignmentFile(f"{self.output_path}/mapped-{i}.bam", "wb", template=in_file, threads=24) as out_file:
                    for record in in_file:
                        record.query_sequence = record.query_name.split("-")[-3]
                        record.template_length = len(record.query_sequence)
                        record.cigartuples = [(0, len(record.query_sequence))]
                        trim_nu = record.query_sequence[-int(record.query_name.split("--")[1]):]
                        record.set_tag("TM", trim_nu)
                        out_file.write(record)
            pysam.sort("-@", "24", "-o", f"{self.output_path}/mapped-{i}.sorted.bam", f"{self.output_path}/mapped-{i}.bam")
            pysam.index(f"{self.output_path}/mapped-{i}.sorted.bam")
            os.remove(f"{self.output_path}/mapped-{i}.sam")
        print(f"[{datetime.datetime.now()}] Merging bam files")
        bam_files = glob.glob(f"{self.output_path}/mapped-*.sorted.bam")
        if len(bam_files) > 0:
            pysam.merge("-@", "24", "-f", f"{self.output_path}/merged-alignment-sorted.sam", *bam_files)
        else:
            print("No bam files found for merging")
    
    def data_process(self):
        print(f"[{datetime.datetime.now()}] Start to remapping the reads!")
        if not os.path.exists(self.output_path):
            os.mkdir(self.output_path)
        self.reformatting_fasta()
        self.remapping()
        self.merge_result()
        print(f"[{datetime.datetime.now()}] Remapping finished!")


##### Argument parsing / help message /version

parser = argparse.ArgumentParser(prog=os.path.basename(__file__))
parser.add_argument('-i', '--input', help='inputdir containing FASTQ files', required=True)
parser.add_argument('-o', '--output', help='output path', required=True)
parser.add_argument('-a', '--adapter', help='whether to trim the adapter and low quality reads', default= "AGATCGGAAGAG")
parser.add_argument('--mir-hairpin', help='reference miRNA hairpin index', 
                    default="/bios-store1/chenyc/Reference_Source/Arabidopsis_Reference/ath_hairpin_bowtie_index/hairpin")
parser.add_argument('--trsno', help='reference trsnoRNA index',
                    default="/bios-store1/chenyc/Reference_Source/Arabidopsis_Reference/ath_trsnoRNA_bowtie_index/ensembl39_araport2016_trsRNA")
parser.add_argument('-v', '--version', action='version', version='%(prog)s 20240307')

args = parser.parse_args()

##### Script control .. go through each file and process

def process_pipeline(input_dir, output_dir, adapter, reference, trsno, flag=True):
    """
    处理整个预处理和比对流程。

    :param input_dir: 输入目录，包含 .fastq.gz 文件
    :param output_dir: 输出目录
    :param adapter: 适配序列
    :param reference: miRNA 的参考序列文件路径
    :param trsno: tRNA/snoRNA 的参考序列文件路径
    :param flag: 临时标志位，默认值为 True
    """
    # 获取所有输入文件
    input_files = glob.glob(os.path.join(input_dir, "*.fastq.gz"))
    if not input_files:
        raise FileNotFoundError(f"No .fastq.gz files found in {input_dir}")

    # 创建必要的输出目录
    format_path = os.path.join(output_dir, "2_formmat_fasta")
    mapping_path = os.path.join(output_dir, "3_remapping")
    os.makedirs(format_path, exist_ok=True)
    os.makedirs(mapping_path, exist_ok=True)

    # 逐个处理文件
    for input_file in input_files:
        sample_name = os.path.basename(input_file).split(".")[0].replace("-", "_")

        # 预处理
        preprocess = PreProcess(input_file, format_path, adapter, sample_name, flag)
        preprocess.data_process()
        print(f"[{datetime.datetime.now()}] Finished preprocessing of {input_file}!")

        # 比对处理
        remapping = ReMapping(mapping_path, reference, trsno, sample_name)
        remapping.data_process()
        print(f"[{datetime.datetime.now()}] Finished remapping of {sample_name}!")

if __name__ == "__main__":
    # # RUN PATH
    # INPUDIR = args.input
    # OUTDIR = args.output
    # FLAG = True # temporary variable, will be removed later
    # INPUT_FILES = glob.glob(os.path.join(INPUDIR, "*.fastq.gz"))

    # FORMAT_PATH = os.path.join(OUTDIR, "2_formmat_fasta")  
    # MAPPING_PATH = os.path.join(OUTDIR, "3_remapping")
    # if not os.path.exists(FORMAT_PATH):
    #     os.mkdir(FORMAT_PATH)
    # if not os.path.exists(MAPPING_PATH):
    #     os.mkdir(MAPPING_PATH)
    # # REF PATH
    # REFERENCE = args.mir_hairpin
    # TRSNORNA = args.trsno
    # ADAPTER = args.adapter

    # for input_file in INPUT_FILES:
    #     SAMPLE_NAME = os.path.basename(input_file).split(".")[0].replace("-", "_")
    #     tag = input_file.split(".")[-1]
        
    #     preprocess = PreProcess(input_file, FORMAT_PATH, ADAPTER, SAMPLE_NAME, FLAG)
    #     preprocess.data_process()
    #     print(f"[{datetime.datetime.now()}] Finish the preprocessing of {input_file}!")

    #     remapping = ReMapping(MAPPING_PATH, REFERENCE, TRSNORNA, SAMPLE_NAME)
    #     remapping.data_process()
    #     print(f"[{datetime.datetime.now()}] Finish the remapping of {SAMPLE_NAME}!")
    process_pipeline(args.input, args.output, args.adapter, args.mir_hairpin, args.trsno)