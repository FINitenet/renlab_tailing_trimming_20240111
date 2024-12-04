#!/usr/bin/python3
# -*- coding:utf-8 -*-
# @FileName  :main.py
# @Time      :2024/11/21 12:00:37
# @Author    :Yuchen@rlab

import logging
from script import (
    tailing_trmming_preprocess,
    tailing_trimming_profile,
    tailing_trimming_length_dist_new,
    tailing_trimming_tail_base_summary_mod,
)

# 设置日志
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


def preprocess_data():
    """步骤1: 数据预处理"""
    logging.info("Starting preprocessing...")
    tailing_trmming_preprocess.preprocess("1_rawdata/", ".")
    logging.info("Preprocessing completed.")


def generate_profile():
    """步骤2: 剪切修饰剖面生成"""
    logging.info("Starting profile generation...")
    tailing_trimming_profile.profile("3_remapping/", ".")
    logging.info("Profile generation completed.")


def run_workflow():
    """步骤3: 运行sRNA测序工作流"""
    logging.info("Starting workflow...")
    workflow_for_srna_seq.workflow("1_rawdata/", ".", 2)
    logging.info("Workflow completed.")


def summarize_mapping_results():
    """步骤4: 汇总比对结果"""
    logging.info("Starting mapping results summary...")
    mapping_results_summary.summarize(".", ".", "4_mapping")
    logging.info("Mapping results summary completed.")


def analyze_length_distribution():
    """步骤5: 分析长度分布"""
    logging.info("Starting length distribution analysis...")
    tailing_trimming_length_dist.analyze("3_remapping/", ".")
    logging.info("Length distribution analysis completed.")


def tail_base_summary():
    """步骤6: Tail base summary"""
    logging.info("Starting tail base summary...")
    tailing_trimming_tail_base_summary_mod.summary(
        filelist_path="/bios-store1/chenyc/Project/YHR/Project_yanghuiru_240510N/filelist",
        project_dir="/bios-store1/chenyc/Project/YHR/Project_yanghuiru_240510N",
        mapping_results="mapping_results_bowtie_4_mapping.csv",
    )
    logging.info("Tail base summary completed.")


def main():
    """主函数"""
    logging.info("Pipeline starting...")
    preprocess_data()
    generate_profile()
    run_workflow()
    summarize_mapping_results()
    analyze_length_distribution()
    tail_base_summary()
    logging.info("Pipeline completed successfully.")


if __name__ == "__main__":
    main()
