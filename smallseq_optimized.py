#!/usr/bin/env python
"""
Optimized SmallSeq Pipeline - Single Script
Processes single-cell small RNA sequencing data (SmallSeq protocol)

Usage:
    python smallseq_pipeline_optimized.py --config config.yaml
    
Or with command line arguments:
    python smallseq_pipeline_optimized.py \
        --rawdata_dir rawdata \
        --output_dir output \
        --genome_dir /path/to/genome \
        --annotation annotations/combined_annots.gp \
        --umi_pattern NNNNNNNN \
        --threads 8

Author: Optimized version
"""

from __future__ import division, print_function
import os
import sys
import argparse
import logging
import subprocess
import pysam
from collections import defaultdict
from multiprocessing import Pool
from functools import partial
import time
import json
import bisect
import gzip

# Try to import optional dependencies
try:
    import yaml
    YAML_AVAILABLE = True
except ImportError:
    YAML_AVAILABLE = False
    print("Warning: PyYAML not available. Config file support disabled.")

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# Hierarchical assignment order from the Small-seq protocol: "Mirbase (miRNAs), GtRNAdb
# (tRNAs), and Gencode transcripts". Lower number wins. Spike-ins and the custom rRNA loci
# sit alongside GENCODE — they never compete with a miRNA or tRNA for the same read.
SOURCE_PRIORITY = {
    'mirbase': 0,
    'gtrnadb': 1,
    'gencode': 2,
    'rrna': 2,
    'spikein': 2,
}


class SmallSeqPipeline:
    """Main pipeline class for SmallSeq data processing"""

    STEP_NAMES = [
        "UMI Removal",
        "Adapter Trimming",
        "STAR Alignment",
        "Soft-clip Removal",
        "Read Length Filtering",
        "UMI Deduplication",
        "Precursor Removal",
        "Count Generation",
        "Long-read Counting",
        "Count Merging",
        "Count Collapsing",
        "Reporting",
    ]

    def __init__(self, config):
        self.config = config
        self.samples = []
        self.checkpoint_file = os.path.join(config['output_dir'], '.pipeline_checkpoint.json')
        self.completed_steps = set()
        self.validate_config()
        self.load_checkpoint()
    
    def load_checkpoint(self):
        """Load completed steps from checkpoint file"""
        if os.path.exists(self.checkpoint_file):
            try:
                with open(self.checkpoint_file, 'r') as f:
                    data = json.load(f)
                    self.completed_steps = set(data.get('completed_steps', []))
                    logger.info(f"Loaded checkpoint with {len(self.completed_steps)} completed steps")
            except Exception as e:
                logger.warning(f"Could not load checkpoint: {e}")
                self.completed_steps = set()
    
    def save_checkpoint(self):
        """Save completed steps to checkpoint file"""
        try:
            self.safe_mkdir(self.config['output_dir'])
            with open(self.checkpoint_file, 'w') as f:
                json.dump({'completed_steps': sorted(list(self.completed_steps))}, f, indent=2)
        except Exception as e:
            logger.warning(f"Could not save checkpoint: {e}")
    
    def mark_step_complete(self, step_name):
        """Mark a step as completed"""
        self.completed_steps.add(step_name)
        self.save_checkpoint()
        logger.info(f"Marked step '{step_name}' as complete")
    
    def is_step_complete(self, step_name):
        """Check if a step has been completed"""
        return step_name in self.completed_steps
    
    def reset_checkpoint(self):
        """Reset all completed steps"""
        self.completed_steps = set()
        if os.path.exists(self.checkpoint_file):
            os.remove(self.checkpoint_file)
        logger.info("Checkpoint reset")
        
    def validate_config(self):
        """Validate configuration parameters"""
        required = ['rawdata_dir', 'output_dir', 'genome_dir', 'annotation','genome_fasta']
        for req in required:
            if req not in self.config:
                raise ValueError(f"Missing required parameter: {req}")
        
        if not os.path.exists(self.config['rawdata_dir']):
            raise FileNotFoundError(f"Raw data directory not found: {self.config['rawdata_dir']}")
        
        # Set defaults
        self.config.setdefault('threads', 4)
        self.config.setdefault('umi_pattern', 'NNNNNNNN')
        
        # Configure read length filtering thresholds
        # 18bp: Minimum from cutadapt adapter trimming (step 2) - NCBI small RNA definition
        # 40bp: Maximum filter (step 6) - Defines upper bound for small RNA molecules
        # 35bp: Precursor threshold (step 8) - Reads ≤35bp always retained (too short to be precursor-derived)
        # 41bp: Reference length minRlen (step 8) - Used for calculating genomic context offset in precursor detection
        self.config.setdefault('max_read_len', 40)
        self.config.setdefault('min_read_len', 41)
        self.config.setdefault('adapter_file', '../adapters/cutadapt_3prime.fa')
        self.config.setdefault('allowed_5p_clip', 0)
        self.config.setdefault('allowed_3p_clip', 3)
        self.config.setdefault('dedup_method', 'adjacency')
        self.config.setdefault('legacy_count', False)
        self.config.setdefault('collapse_level', 'gene')

    def safe_mkdir(self, path):
        """Create directory if it doesn't exist"""
        if not os.path.exists(path):
            os.makedirs(path, mode=0o774)
            logger.info(f"Created directory: {path}")

    def get_untrimmed_max_len(self):
        """Longest read seen after UMI+CA removal but before length filtering.

        The protocol derives its length thresholds from the sequencing read length: with a
        51bp run, 51 - 8 (UMI) - 2 (CA) = 41, hence maxRlen=40/minRlen=41. Step 8's precursor
        check reconstructs "the read was this long before adapter trimming", so it needs the
        REAL value, not min_read_len — which is only correct for a 51bp run. Detect it once
        from the step-2 FASTQs and cache it.
        """
        if self.config.get('untrimmed_max_len'):
            return self.config['untrimmed_max_len']

        trimmed_dir = os.path.join(self.config['output_dir'], 'step2_adapter_trimmed')
        maxlen = 0
        for sample in self.samples:
            fq = os.path.join(trimmed_dir, sample, f"{sample}.fastq.gz")
            if not os.path.exists(fq):
                continue
            with gzip.open(fq, 'rt') as fh:
                for i, line in enumerate(fh):
                    if i % 4 == 1:
                        maxlen = max(maxlen, len(line.strip()))
                    if i > 400000:  # a couple of hundred thousand reads is plenty
                        break
            if maxlen:
                break

        if not maxlen:
            maxlen = self.config['min_read_len']
            logger.warning(f"Could not detect untrimmed read length; falling back to "
                           f"min_read_len={maxlen}. Step 8 may filter the wrong lengths.")
        else:
            logger.info(f"Detected untrimmed max read length: {maxlen}nt "
                        f"(sequencing length minus 8nt UMI and 2nt CA)")
            if maxlen != self.config['min_read_len']:
                logger.warning(
                    f"min_read_len={self.config['min_read_len']} but the real untrimmed length "
                    f"is {maxlen}nt. Step 8's precursor check will use {maxlen}; reads longer "
                    f"than max_read_len={self.config['max_read_len']} go to step6_long/.")

        self.config['untrimmed_max_len'] = maxlen
        return maxlen


    def run(self):
        """Execute the complete pipeline"""
        logger.info("="*60)
        logger.info("Starting SmallSeq Pipeline")
        logger.info("="*60)
        
        start_time = time.time()
        
        # Get sample list
        self.samples = [s for s in os.listdir(self.config['rawdata_dir']) 
                       if os.path.isdir(os.path.join(self.config['rawdata_dir'], s))]
        logger.info(f"Found {len(self.samples)} samples to process")
        
        # Pipeline steps
        step_funcs = [
            self.step1_remove_umi,
            self.step2_trim_adapters,
            self.step3_star_alignment,
            self.step5_remove_softclipped,
            self.step6_filter_by_length,
            self.step7_umi_dedup,
            self.step8_remove_precursors,
            self.step9_count_smallrnas,
            self.step10_count_long,
            self.step11_merge_counts,
            self.step12_collapse_counts,
            self.step13_reporting,
        ]
        steps = list(zip(self.STEP_NAMES, step_funcs))
        
        for step_name, step_func in steps:
            if self.is_step_complete(step_name):
                logger.info(f"{'='*60}")
                logger.info(f"Step: {step_name} [SKIPPED - Already completed]")
                logger.info(f"{'='*60}")
                continue
            
            logger.info(f"\n{'='*60}")
            logger.info(f"Step: {step_name}")
            logger.info(f"{'='*60}")
            try:
                step_func()
                self.mark_step_complete(step_name)
            except Exception as e:
                logger.error(f"Error in {step_name}: {str(e)}")
                raise
        
        elapsed = time.time() - start_time
        logger.info(f"\n{'='*60}")
        logger.info(f"Pipeline completed successfully in {elapsed:.2f} seconds")
        logger.info(f"Output: {os.path.join(self.config['output_dir'], 'counts_molc_final.txt')}")
        logger.info(f"{'='*60}")
    
    # ===== Step 1: UMI Removal =====
    def _process_umi_sample(self, sample):
        """Process a single sample for UMI removal"""
        input_dir = self.config['rawdata_dir']
        output_dir = os.path.join(self.config['output_dir'], 'step1_umi_removed')
        
        sample_in = os.path.join(input_dir, sample)
        sample_out = os.path.join(output_dir, sample)
        self.safe_mkdir(sample_out)
        
        # Find FASTQ file
        fq_files = [f for f in os.listdir(sample_in) if f.endswith(('.fq', '.fastq', '.fq.gz', '.fastq.gz'))]
        if not fq_files:
            logger.warning(f"No FASTQ file found for {sample}")
            return
        
        raw_fq = os.path.join(sample_in, fq_files[0])
        fq_name = fq_files[0].split('.')[0]
        trimmed_fq = os.path.join(sample_out, f"{fq_name}_umiTrim.fq.gz")
        logfile = os.path.join(sample_out, "extract.log")
        
        # Use UMI-tools
        cmd = f"umi_tools extract --bc-pattern={self.config['umi_pattern']} " \
              f"-L {logfile} -I {raw_fq} -S {trimmed_fq}"
        
        result = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
        if result.returncode != 0:
            logger.warning(f"UMI extraction warning for {sample}: {result.stderr}")
        #result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        #if result.returncode != 0:
        #    logger.warning(f"UMI extraction warning for {sample}: {result.stderr}")
    
    def step1_remove_umi(self):
        """Remove UMI sequences from FASTQ files"""
        output_dir = os.path.join(self.config['output_dir'], 'step1_umi_removed')
        self.safe_mkdir(output_dir)
        
        with Pool(self.config['threads']) as pool:
            pool.map(self._process_umi_sample, self.samples)
    
    # ===== Step 2: Adapter Trimming =====
    def _process_adapter_sample(self, sample):
        """Process a single sample for adapter trimming"""
        input_dir = os.path.join(self.config['output_dir'], 'step1_umi_removed')
        output_dir = os.path.join(self.config['output_dir'], 'step2_adapter_trimmed')
        
        sample_out = os.path.join(output_dir, sample)
        self.safe_mkdir(sample_out)
        
        input_fq = os.path.join(input_dir, sample, f"{sample}_umiTrim.fq.gz")
        output_fq = os.path.join(sample_out, f"{sample}.fastq.gz")
        logfile = os.path.join(sample_out, "cutadapt.log")

        if not os.path.exists(input_fq):
            logger.warning(f"Input file not found: {input_fq}")
            return

        cmd = f"cutadapt -a file:{self.config['adapter_file']} " \
              f"-e 0.1 -O 1 -u 2 --minimum-length 18 " \
              f"-o {output_fq} {input_fq}"

        #subprocess.run(cmd, shell=True, check=True)
        result = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
        with open(logfile, 'w') as fh:
            fh.write(result.stdout)
        if result.returncode != 0:
            logger.warning(f"Adapter trimming warning for {sample}: {result.stderr}")
    
    def step2_trim_adapters(self):
        """Trim adapters using cutadapt"""
        output_dir = os.path.join(self.config['output_dir'], 'step2_adapter_trimmed')
        self.safe_mkdir(output_dir)
        
        with Pool(self.config['threads']) as pool:
            pool.map(self._process_adapter_sample, self.samples)
    
    # ===== Step 3: STAR Alignment =====
    def step3_star_alignment(self):
        """Run STAR once with read groups and split BAMs per sample"""
        output_dir = os.path.join(self.config['output_dir'], 'step3_star_aligned')
        input_dir = os.path.join(self.config['output_dir'], 'step2_adapter_trimmed')
        self.safe_mkdir(output_dir)
        
        fq_list = os.path.abspath(os.path.join(output_dir, "star_fastqs.txt"))
        rg_list = os.path.abspath(os.path.join(output_dir, "star_readgroups.txt"))
        samples_used = []
        
        # 1) Collect FASTQs + read groups
        with open(fq_list, "w") as fq_fh, open(rg_list, "w") as rg_fh:
            for sample in self.samples:
                fq = os.path.abspath(os.path.join(input_dir, sample, f"{sample}.fastq.gz"))
                if not os.path.exists(fq):
                    logger.warning(f"Input file not found: {fq}")
                    continue
                fq_fh.write(f"{fq}\t-\tID:{sample}\tSM:{sample}\tPL:ILLUMINA\n")
                rg_fh.write(f"ID:{sample} SM:{sample} PL:ILLUMINA\n")
                samples_used.append(sample)
        
        if not samples_used:
            raise RuntimeError("No valid FASTQs found for STAR alignment")
        
        # 2) Run STAR
        genome_dir = os.path.abspath(self.config['genome_dir'])
        prefix = os.path.abspath(os.path.join(output_dir, "all_samples_"))
        star_cmd = (
            f"STAR "
            f"--runThreadN {self.config['threads']} "
            f"--genomeDir {genome_dir} "
            f"--readFilesManifest {fq_list} "
            f"--readFilesCommand zcat "
            f"--outSAMtype BAM Unsorted "
            f"--outSAMstrandField intronMotif "
            f"--outFilterMultimapNmax 20 "
            f"--outFilterScoreMinOverLread 0 "
            f"--outFilterMatchNmin 18 "
            f"--outFilterMatchNminOverLread 0 "
            f"--outFilterMismatchNoverLmax 0.04 "
            f"--alignIntronMax 1 "
            f"--outFileNamePrefix {prefix} "
            f"--outSAMattributes NH HI AS nM RG"
        )
        logger.info(f"Running STAR aligmnet : {star_cmd}")
        result = subprocess.run(
            star_cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        if result.returncode != 0:
            raise RuntimeError(f"STAR failed:\n{result.stderr}")
        
        # 3a) Sort the unsorted STAR output, then remove it immediately to free space
        input_combined_bam = os.path.abspath(os.path.join(output_dir, "all_samples_Aligned.out.bam"))
        combined_bam = os.path.abspath(os.path.join(output_dir, "all_samples_Aligned.sortedByCoord.out.bam"))

        if not os.path.exists(input_combined_bam):
            raise FileNotFoundError("Combined STAR BAM not found")

        sort_cmd = f"samtools sort -@ {self.config['threads']} {input_combined_bam} -o {combined_bam}"
        logger.info(f"Sorting combined BAM: {sort_cmd}")
        result = subprocess.run(
            sort_cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        if result.returncode != 0:
            raise RuntimeError(f"samtools sort failed:\n{result.stderr}")

        os.remove(input_combined_bam)
        logger.info(f"Removed unsorted BAM: {input_combined_bam}")

        # 3b) Split sorted BAM by read group, then remove it to free space
        split_cmd = (
            f"samtools split -@ {self.config['threads']} -M -1 -d RG "
            f"-f {os.path.abspath(output_dir)}/%*_%!.%. {combined_bam}"
        )
        logger.info(f"Splitting sorted BAM: {split_cmd}")
        result = subprocess.run(
            split_cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        if result.returncode != 0:
            raise RuntimeError(f"samtools split failed:\n{result.stderr}")

        os.remove(combined_bam)
        logger.info(f"Removed sorted combined BAM: {combined_bam}")
        
        # 4) Move, rename, and index per sample
        for sample in samples_used:
            rg_bam = os.path.join(
                output_dir,
                f"all_samples_Aligned.sortedByCoord.out_{sample}.bam",
            )
            if not os.path.exists(rg_bam):
                logger.warning(f"No BAM produced for sample {sample}")
                continue
            
            sample_dir = os.path.join(output_dir, sample)
            self.safe_mkdir(sample_dir)
            final_bam = os.path.join(sample_dir, f"{sample}.bam")
            os.rename(rg_bam, final_bam)
            
            subprocess.run(
                f"samtools index {final_bam}",
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
            )
    
    # ===== Step 4: SAM Processing =====
    #def _process_sam_sample(self, sample):
    #    """Process a single sample SAM file"""
    #    input_dir = os.path.join(self.config['output_dir'], 'step3_star_aligned')
    #    
    #    sample_dir = os.path.join(input_dir, sample)
    #    sam_file = os.path.join(sample_dir, "Aligned.out.sam")
    #    bam_file = os.path.join(sample_dir, f"{sample}.bam")
    #    
    #    if not os.path.exists(sam_file):
    #        logger.warning(f"SAM file not found: {sam_file}")
    #        return
    #    
    #    # Convert to BAM and sort
    #    pysam.view('-bS', sam_file, '-o', bam_file + '.tmp', catch_stdout=False)
    #    pysam.sort('-o', bam_file, bam_file + '.tmp')
    #    pysam.index(bam_file)
    #    
    #    # Clean up
    #    os.remove(bam_file + '.tmp')
    #    os.remove(sam_file)  # Save space
    #
    #def step4_process_sam(self):
    #    """Convert SAM to sorted BAM and index"""
    #    with Pool(self.config['threads']) as pool:
    #        pool.map(self._process_sam_sample, self.samples)
    
    # ===== Step 5: Soft-clip Removal =====
    def _process_softclip_sample(self, sample):
        """Process a single sample for soft-clip removal"""
        input_dir = os.path.join(self.config['output_dir'], 'step3_star_aligned')
        output_dir = os.path.join(self.config['output_dir'], 'step5_clipped_removed')
        
        sample_out = os.path.join(output_dir, sample)
        self.safe_mkdir(sample_out)
        
        inbam = os.path.join(input_dir, sample, f"{sample}.bam")
        outbam_tmp = os.path.join(sample_out, f"{sample}_tmp.bam")
        outbam = os.path.join(sample_out, f"{sample}.bam")
        
        if not os.path.exists(inbam):
            logger.warning(f"Input BAM not found: {inbam}")
            return
        
        inbam_obj = pysam.AlignmentFile(inbam, "rb")
        outbam_obj = pysam.AlignmentFile(outbam_tmp, "wb", template=inbam_obj)
        
        total, kept = 0, 0
        for read in inbam_obj:
            total += 1
            cigar = read.cigar
            
            # Check clipping
            clip_5p = clip_3p = hardclip = ins = dels = 0
            for i, (op, length) in enumerate(cigar):
                if op == 5: hardclip += length
                elif op == 1: ins += length
                elif op == 2: dels += length
                elif i == 0 and op == 4: clip_5p += length
                elif i == len(cigar)-1 and op == 4: clip_3p += length
            
            # Filter
            if (clip_5p <= self.config['allowed_5p_clip'] and 
                clip_3p <= self.config['allowed_3p_clip'] and
                hardclip == 0 and ins == 0 and dels == 0):
                outbam_obj.write(read)
                kept += 1
        
        outbam_obj.close()
        inbam_obj.close()
        
        # Sort and index
        pysam.sort('-o', outbam, outbam_tmp)
        pysam.index(outbam)
        os.remove(outbam_tmp)

        logger.info(f"{sample}: Removed {100*(1-kept/total):.2f}% clipped reads")

        stats_file = os.path.join(sample_out, f"{sample}_softclip_stats.txt")
        with open(stats_file, 'w') as fh:
            fh.write(f"total_reads\t{total}\n")
            fh.write(f"kept_reads\t{kept}\n")
            fh.write(f"removed_pct\t{100*(1-kept/total):.2f}\n")
    
    def step5_remove_softclipped(self):
        """Remove soft-clipped reads"""
        output_dir = os.path.join(self.config['output_dir'], 'step5_clipped_removed')
        self.safe_mkdir(output_dir)
        
        with Pool(self.config['threads']) as pool:
            pool.map(self._process_softclip_sample, self.samples)
    
    # ===== Step 6: Read Length Filtering =====
    def _process_readlen_sample(self, sample):
        """Process a single sample for read length filtering"""
        input_dir = os.path.join(self.config['output_dir'], 'step5_clipped_removed')
        output_dir = os.path.join(self.config['output_dir'], f'step6_max{self.config["max_read_len"]}')
        
        sample_out = os.path.join(output_dir, sample)
        self.safe_mkdir(sample_out)
        
        long_dir = os.path.join(self.config['output_dir'], 'step6_long')
        long_out = os.path.join(long_dir, sample)
        self.safe_mkdir(long_out)

        inbam = os.path.join(input_dir, sample, f"{sample}.bam")
        outbam = os.path.join(sample_out, f"{sample}.bam")
        longbam = os.path.join(long_out, f"{sample}.bam")

        maxlen = self.config['max_read_len']

        # Reads longer than the cap are NOT small RNAs by the protocol's definition, but they
        # are real data (tRNA/snoRNA fragments) — route them to step6_long/ instead of
        # discarding. With 101bp sequencing this is >50% of the library.
        cmd = f"samtools view -h {inbam} | " \
              f"awk 'length($10) <= {maxlen} || $1~\"@\"' | " \
              f"samtools view -bS - > {outbam}"
        result = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
        if result.returncode !=0:
            logger.warning(f"Filtering warning for {sample}: {result.stderr}")

        long_cmd = f"samtools view -h {inbam} | " \
                   f"awk 'length($10) > {maxlen} || $1~\"@\"' | " \
                   f"samtools view -bS - > {longbam}"
        result = subprocess.run(long_cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
        if result.returncode != 0:
            logger.warning(f"Long-read split warning for {sample}: {result.stderr}")

        pysam.index(outbam)
        pysam.index(longbam)

        n_short = int(pysam.view("-c", outbam).strip() or 0)
        n_long = int(pysam.view("-c", longbam).strip() or 0)
        logger.info(f"{sample}: {n_short} alignments <={maxlen}nt -> step6_max{maxlen}, "
                    f"{n_long} >{maxlen}nt -> step6_long")

        input_flagstat = os.path.join(sample_out, f"{sample}_input_flagstat.txt")
        output_flagstat = os.path.join(sample_out, f"{sample}_output_flagstat.txt")
        for bam, flagstat_file in ((inbam, input_flagstat), (outbam, output_flagstat)):
            result = subprocess.run(f"samtools flagstat {bam} > {flagstat_file}",
                                    shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
            if result.returncode != 0:
                logger.warning(f"flagstat warning for {sample}: {result.stderr}")
    
    def step6_filter_by_length(self):
        """Split reads by length: <=max_read_len are small RNAs, longer ones are kept aside

        Long reads are not discarded — they go to step6_long/ for the separate tRNA/MINTmap
        pipeline. Nothing downstream in this pipeline consumes them.
        """
        output_dir = os.path.join(self.config['output_dir'], f'step6_max{self.config["max_read_len"]}')
        self.safe_mkdir(output_dir)
        self.safe_mkdir(os.path.join(self.config['output_dir'], 'step6_long'))

        with Pool(self.config['threads']) as pool:
            pool.map(self._process_readlen_sample, self.samples)
    
    # ===== Step 7: UMI Deduplication =====
    def _process_dedup_sample(self, sample, stage=''):
        """Process a single sample for UMI deduplication

        stage='' processes the small-RNA path (step6_max{N} -> step7_dedup).
        stage='_long' processes the >max_read_len path (step6_long -> step7_dedup_long).
        """
        if stage == '_long':
            input_dir = os.path.join(self.config['output_dir'], 'step6_long')
        else:
            input_dir = os.path.join(self.config['output_dir'], f'step6_max{self.config["max_read_len"]}')
        output_dir = os.path.join(self.config['output_dir'], f'step7_dedup{stage}')

        sample_out = os.path.join(output_dir, sample)
        self.safe_mkdir(sample_out)
        
        inbam = os.path.join(input_dir, sample, f"{sample}.bam")
        outbam = os.path.join(sample_out, f"{sample}_dedup.bam")
        logfile = os.path.join(sample_out, "dedup.log")
        
        # NOTE: --read-length flag is CRITICAL for small RNA-seq. It ensures reads of different
        # lengths at the same position with the same UMI are kept as separate molecules.
        # Without it, miR-21 (22bp) and miR-21-5p (23bp) would be incorrectly merged.
        cmd = f"umi_tools dedup --method {self.config['dedup_method']} " \
              f"--read-length " \
              f"--output-stats {sample_out}/stats " \
              f"-I {inbam} -S {outbam} -L {logfile}"
        
        result = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
        if result.returncode !=0:
            logger.warning(f"UMI deduplication warning for {sample}: {result.stderr}")
    
    def step7_umi_dedup(self):
        """Remove PCR duplicates using UMI"""
        output_dir = os.path.join(self.config['output_dir'], 'step7_dedup')
        self.safe_mkdir(output_dir)

        with Pool(self.config['threads']) as pool:
            pool.map(self._process_dedup_sample, self.samples)
    
    # ===== Step 8: Precursor Removal =====
    def _process_precursor_sample(self, sample, gf, stage=''):
        """Process a single sample for precursor removal

        stage='' is the small-RNA path; stage='_long' is the >max_read_len path. The check is
        most meaningful on the long path: those are the max-length reads it was designed for.
        """
        input_dir = os.path.join(self.config['output_dir'], f'step7_dedup{stage}')
        output_dir = os.path.join(self.config['output_dir'], f'step8_precursor_removed{stage}')

        sample_out = os.path.join(output_dir, sample)
        self.safe_mkdir(sample_out)
        
        inbam = os.path.join(input_dir, sample, f"{sample}_dedup.bam")
        outbam_tmp = os.path.join(sample_out, f"{sample}_tmp.bam")
        outbam = os.path.join(sample_out, f"{sample}.bam")
        
        inbam_obj = pysam.AlignmentFile(inbam, "rb")
        outbam_obj = pysam.AlignmentFile(outbam_tmp, "wb", template=inbam_obj)

        # A read only needs the precursor check if adapter trimming could plausibly have
        # eaten just a base or two off a FULL-LENGTH read. That means keying on the real
        # untrimmed read length, not min_read_len (which is only equal to it for a 51bp run).
        # Using min_read_len on longer sequencing runs makes this filter fire on genuine
        # small RNAs and destroy them.
        untrimmed = self.config['untrimmed_max_len']

        total, kept = 0, 0
        for read in inbam_obj:
            total += 1
            readlen = len(read.query_sequence)

            if readlen <= 35:  # Keep short reads
                outbam_obj.write(read)
                kept += 1
                continue

            readchr = inbam_obj.get_reference_name(read.reference_id)
            readstart = read.pos + 1
            readend = read.reference_end

            upperlimit = untrimmed - readlen

            # Keep by default. Only an actual adapter-lookalike in the genome removes a read.
            # (The original keyed on flag == 0 / 16, which silently dropped every secondary
            # alignment of a read longer than 35nt and broke multi-mapper weighting.)
            keep_read = True
            if 0 < upperlimit <= 5:
                if not read.is_reverse:
                    bpwindow = gf.get_seq_from_to(readchr, readend+1, readend+upperlimit)
                    patterns = {
                        untrimmed-1: ["T", "A"],
                        untrimmed-2: ["TG", "AA"],
                        untrimmed-3: ["TGG", "AAA"],
                        untrimmed-4: ["TGGA", "AAAA"],
                        untrimmed-5: ["TGGAA", "AAAAA"]
                    }
                else:
                    bpwindow = gf.get_seq_from_to(readchr, readstart-upperlimit, readstart-1)
                    patterns = {
                        untrimmed-1: ["A", "T"],
                        untrimmed-2: ["CA", "TT"],
                        untrimmed-3: ["CCA", "TTT"],
                        untrimmed-4: ["TCCA", "TTTT"],
                        untrimmed-5: ["TTCCA", "TTTTT"]
                    }
                if readlen in patterns and bpwindow in patterns[readlen]:
                    keep_read = False

            if keep_read:
                outbam_obj.write(read)
                kept += 1

        outbam_obj.close()
        inbam_obj.close()

        # Sort and index
        pysam.sort('-o', outbam, outbam_tmp)
        pysam.index(outbam)
        os.remove(outbam_tmp)

        stats_file = os.path.join(sample_out, f"{sample}_precursor_stats.txt")
        with open(stats_file, 'w') as fh:
            fh.write(f"total_reads\t{total}\n")
            fh.write(f"kept_reads\t{kept}\n")
            fh.write(f"removed_pct\t{100*(1-kept/total):.2f}\n")

    def _process_precursor_sample_parallel(self, sample, stage=''):
        """Worker-safe wrapper: creates its own GenomeFetch instance per process"""
        sys.path.insert(0, os.path.dirname(__file__))
        from GenomeFetch import GenomeFetch
        gf = GenomeFetch(genomedir=self.config['genome_fasta'])
        self._process_precursor_sample(sample, gf, stage=stage)

    def step8_remove_precursors(self):
        """Remove reads from precursor RNAs based on genomic context"""
        output_dir = os.path.join(self.config['output_dir'], 'step8_precursor_removed')
        self.safe_mkdir(output_dir)

        # Resolve before forking so every worker inherits the same value
        self.get_untrimmed_max_len()

        with Pool(self.config['threads']) as pool:
            pool.map(self._process_precursor_sample_parallel, self.samples)
    
    # ===== Step 9: Count Small RNAs =====
    def _process_count_sample(self, sample, interval_lists, interval_starts, max_interval_len,
                              coord2geneid, geneid2name, geneidlist, legacy_count,
                              stage='', exclude_prios=frozenset()):
        """Process a single sample for counting using gene annotations

        This function is designed to be called safely in parallel via multiprocessing.Pool.
        It uses a locally-scoped midpos variable, making it thread-safe (unlike the original
        count_smallrnas.py which had a global midpos variable bug that forced serial execution).

        exclude_prios: source priorities that may not be assigned. The long path passes
        miRBase's priority: mature miRNAs are 16-28nt, so a 41-91nt read can never BE one, but
        the overlap test is an intersection and would happily match a read that merely grazes
        one -- and miRBase outranks GENCODE, so that graze would discard the correct host-gene
        assignment.
        """
        input_dir = os.path.join(self.config['output_dir'], f'step8_precursor_removed{stage}')
        output_dir = os.path.join(self.config['output_dir'], f'step9_counts{stage}')

        sample_out = os.path.join(output_dir, sample)
        self.safe_mkdir(sample_out)

        inbam = os.path.join(input_dir, sample, f"{sample}.bam")
        outfile = os.path.join(sample_out, f"{sample}_Count.txt")

        def find_overlaps(chrom, qstart, qend, strand):
            """Find gene intervals overlapping [qstart, qend] using bisect O(log n) lookup.
            For legacy midpoint mode, pass qstart == qend == midpos.

            Applies the protocol's hierarchical assignment: when a read overlaps annotations
            from several databases, only the highest-priority source is kept
            (miRBase > GtRNAdb > GENCODE). Ties within that source are returned together so
            the caller's 1/annotatedCount weighting still applies.
            """
            key = (chrom, strand)
            if key not in interval_lists:
                return []
            ivs = interval_lists[key]
            starts = interval_starts[key]
            i = bisect.bisect_right(starts, qend) - 1
            hits = []
            best = None
            while i >= 0 and starts[i] >= qstart - max_interval_len:
                start, end, geneid, prio = ivs[i]
                # Skipped before the priority reduction, so an excluded source can neither be
                # assigned nor suppress a legitimate lower-priority hit.
                if end >= qstart and prio not in exclude_prios:
                    hits.append((prio, f"{chrom}:{start+1}-{end}:{strand}"))
                    if best is None or prio < best:
                        best = prio
                i -= 1
            if best is None:
                return []
            return [coord for prio, coord in hits if prio == best]

        inbam_obj = pysam.AlignmentFile(inbam, "rb")

        read2overlaps = defaultdict(list)

        for read in inbam_obj:
            readchr = inbam_obj.get_reference_name(read.reference_id)
            readstart = read.pos
            readend = read.reference_end
            strand = "-" if read.is_reverse else "+"

            if legacy_count:
                midpos = (readstart + readend) // 2
                overlaps = find_overlaps(readchr, midpos, midpos, strand)
            else:
                overlaps = find_overlaps(readchr, readstart, readend, strand)
            read2overlaps[read.qname].append(overlaps)
        
        inbam_obj.close()
        
        # Count
        geneid2counts = {}
        num_unannot = 0
        
        for read, overlap_list in read2overlaps.items():
            read_count = len(overlap_list)
            annot_count = sum(1 for ol in overlap_list if ol)
            
            if annot_count > 0:
                for overlaps in overlap_list:
                    if overlaps:
                        # Take only the first overlapping coordinate,
                        # matching original count_smallrnas.py coord[0] behavior
                        coord = overlaps[0]
                        geneid = coord2geneid.get(coord, 'NA')
                        if geneid not in geneid2counts:
                            geneid2counts[geneid] = 0
                        geneid2counts[geneid] += 1 / annot_count
                    else:
                        # Unannotated alignment of a partially-annotated read
                        # still contributes 1/annot_count to "NA", matching original
                        if 'NA' not in geneid2counts:
                            geneid2counts['NA'] = 0
                        geneid2counts['NA'] += 1 / annot_count
            else:
                num_unannot += 1
        
        num_annot = sum(v for k, v in geneid2counts.items() if k != 'NA' and not k.startswith('P-cel'))
        
        # Write output. The long pass labels its column {sample}_long so it stays distinct
        # from the small-RNA column when the two are merged side by side.
        col_name = f"{sample}_long" if stage == '_long' else sample
        with open(outfile, 'w') as fh:
            fh.write(f"#samples\t{col_name}\n")
            fh.write(f"#unannotatedmolc\t{num_unannot}\n")
            fh.write(f"#annotatedmolc\t{num_annot}\n")
            for geneid in geneidlist:
                fh.write(f"{geneid2name[geneid]}\t{geneid}\t{geneid2counts.get(geneid, 0)}\n")
    
    def _load_annotation(self):
        """Parse the GenePred annotation into bisect-ready per-(chrom,strand) exon intervals.

        Shared by the small-RNA and long-read counting passes.
        """
        geneid2name = {}
        coord2geneid = {}
        geneidlist = []

        # interval_lists: (chrom, strand) -> sorted list of (start, end, geneid)
        # interval_starts: (chrom, strand) -> sorted list of starts (parallel, for bisect)
        interval_lists = {}
        max_interval_len = 0

        missing_source = False
        for line in open(self.config['annotation'], 'r'):
            p = line.split()
            chrom, strand, geneid, genename = p[2], p[3], p[1], p[12]
            # Col 16 is the source database, written by build_annotation.py. Older .gp files
            # do not have it; fall back to a flat priority so they still run.
            if len(p) > 16:
                source = p[16]
            else:
                source = 'gencode'
                missing_source = True
            prio = SOURCE_PRIORITY.get(source, SOURCE_PRIORITY['gencode'])

            key = (chrom, strand)
            if key not in interval_lists:
                interval_lists[key] = []

            # One interval per EXON, not one spanning txStart..txEnd. A transcript-span
            # window makes every intron part of the gene's counting window, which is how
            # intronic reads were being credited to protein-coding genes.
            exon_starts = [int(x) for x in p[9].rstrip(',').split(',')]
            exon_ends = [int(x) for x in p[10].rstrip(',').split(',')]
            for start, end in zip(exon_starts, exon_ends):
                interval_lists[key].append((start, end, geneid, prio))
                max_interval_len = max(max_interval_len, end - start)
                coord = f"{chrom}:{start+1}-{end}:{strand}"
                coord2geneid[coord] = geneid

            geneid2name[geneid] = genename
            geneidlist.append(geneid)

        if missing_source:
            logger.warning(
                "Annotation has no source column (col 17) — hierarchical assignment "
                "(miRBase > GtRNAdb > GENCODE) is DISABLED. Regenerate with build_annotation.py.")

        # Sort by start and build parallel start lists for bisect
        interval_starts = {}
        for key in interval_lists:
            interval_lists[key].sort()
            interval_starts[key] = [iv[0] for iv in interval_lists[key]]

        return dict(interval_lists=interval_lists,
                    interval_starts=interval_starts,
                    max_interval_len=max_interval_len,
                    coord2geneid=coord2geneid,
                    geneid2name=geneid2name,
                    geneidlist=geneidlist)

    def _run_counting(self, stage='', exclude_prios=frozenset()):
        """Count one BAM stage against the annotation, in parallel across samples."""
        self.safe_mkdir(os.path.join(self.config['output_dir'], f'step9_counts{stage}'))
        annot = self._load_annotation()
        count_func = partial(self._process_count_sample,
                             legacy_count=self.config['legacy_count'],
                             stage=stage,
                             exclude_prios=exclude_prios,
                             **annot)
        with Pool(self.config['threads']) as pool:
            pool.map(count_func, self.samples)

    def step9_count_smallrnas(self):
        """Count small RNAs using gene annotation file

        PARALLELIZATION NOTE: The original pipeline's count_smallrnas.py had parallelization
        disabled due to a global 'midpos' variable bug that corrupted results in parallel
        (see src/count_smallrnas.py line 139 comment). This version fixes that bug by using
        locally-scoped variables in _process_count_sample(), enabling safe parallel execution
        via multiprocessing.Pool. This provides significant performance improvement over the
        original serial implementation.
        """
        self._run_counting()

    # ===== Step 10: Count the long-read fraction =====
    def step10_count_long(self):
        """Dedup, precursor-filter and count the >max_read_len reads set aside by step 6.

        These reads are NOT small RNAs by the protocol's definition, but they are real data
        (tRNA fragments are 55-176nt, snoRNAs longer). They run the same 7/8/9 path as the
        small-RNA fraction so their counts are UMI-deduplicated MOLECULES on the same footing,
        with one difference: miRBase is excluded as an assignment target (see
        _process_count_sample).
        """
        long_dir = os.path.join(self.config['output_dir'], 'step6_long')
        if not os.path.isdir(long_dir):
            logger.info("No step6_long/ directory; skipping long-read counting.")
            return

        self.safe_mkdir(os.path.join(self.config['output_dir'], 'step7_dedup_long'))
        with Pool(self.config['threads']) as pool:
            pool.map(partial(self._process_dedup_sample, stage='_long'), self.samples)

        self.safe_mkdir(os.path.join(self.config['output_dir'], 'step8_precursor_removed_long'))
        self.get_untrimmed_max_len()
        with Pool(self.config['threads']) as pool:
            pool.map(partial(self._process_precursor_sample_parallel, stage='_long'), self.samples)

        self._run_counting(stage='_long',
                           exclude_prios=frozenset({SOURCE_PRIORITY['mirbase']}))


    # ===== Step 11: Merge Counts =====
    def step11_merge_counts(self):
        """Merge count files from all samples

        Emits the small-RNA column per sample followed by, where present, that sample's
        long-read column ({sample}_long). The long columns are purely additive: the small-RNA
        columns are byte-identical to a run without long-read counting.
        """
        input_dir = os.path.join(self.config['output_dir'], 'step9_counts')
        long_dir = os.path.join(self.config['output_dir'], 'step9_counts_long')
        output_file = os.path.join(self.config['output_dir'], 'counts_molc.txt')

        count_files = []
        for sample in self.samples:
            count_file = os.path.join(input_dir, sample, f"{sample}_Count.txt")
            if os.path.exists(count_file):
                count_files.append(count_file)
        for sample in self.samples:
            long_file = os.path.join(long_dir, sample, f"{sample}_Count.txt")
            if os.path.exists(long_file):
                count_files.append(long_file)

        if not count_files:
            logger.error("No count files found!")
            return
        
        # Parse files
        header = ['#samples', '#unannotatedmolc', '#annotatedmolc']
        genelines = []
        
        for i, inf in enumerate(count_files):
            with open(inf, 'r') as fh:
                gene_idx = 0
                for line in fh:
                    p = line.strip().split('\t')
                    if p[0] == '#samples':
                        header[0] += '\t' + '\t'.join(p[1:])
                    elif p[0] == '#unannotatedmolc':
                        header[1] += '\t' + '\t'.join(p[1:])
                    elif p[0] == '#annotatedmolc':
                        header[2] += '\t' + '\t'.join(p[1:])
                    elif not line.startswith('#'):
                        if i == 0:
                            genelines.append('\t'.join(p[:2]))
                        genelines[gene_idx] += '\t' + '\t'.join(p[2:])
                        gene_idx += 1
        
        # Write output
        with open(output_file, 'w') as fh:
            for h in header:
                fh.write(h + '\n')
            for line in genelines:
                fh.write(line + '\n')
        
        logger.info(f"Merged counts written to {output_file}")
    
    # ===== Step 12: Collapse counts =====
    def step12_collapse_counts(self):
        """Sum counts onto gene-level rows (or per-transcript rows in legacy mode)

        Reads are assigned to exactly one transcript ID by step 9, and transcripts of a
        gene share exons, so which transcript of a gene a read lands on is arbitrary.
        Summing rows that share a gene name is the meaningful aggregation.

        collapse_level:
          'gene'       - sum every row sharing a gene name, including across loci
                         (multi-copy families like Y_RNA/U6 are indistinguishable by
                         short reads, so their copies merge into one row)
          'transcript' - legacy behaviour: report one row per transcript, collapsing
                         only miRBase mature miRNAs, which merge across their loci
        """
        input_file = os.path.join(self.config['output_dir'], 'counts_molc.txt')
        output_file = os.path.join(self.config['output_dir'], 'counts_molc_final.txt')
        level = self.config['collapse_level']

        gene2counts = {}   # gene name -> element-wise summed counts
        gene2ntx = {}      # gene name -> number of transcript rows summed into it
        order = []         # first-appearance order, so output is deterministic

        with open(output_file, 'w') as outfh:
            for line in open(input_file, 'r'):
                if line.startswith('#'):
                    outfh.write(line)
                    continue

                p = line.rstrip('\n').split('\t')
                genename = p[0]

                collapse = (level == 'gene') or genename.startswith('hsa')
                if not collapse:
                    outfh.write(line)
                    continue

                counts = [float(c) for c in p[2:]]
                if genename not in gene2counts:
                    gene2counts[genename] = [0.0] * len(counts)
                    gene2ntx[genename] = 0
                    order.append(genename)
                gene2counts[genename] = [a + b for a, b in zip(gene2counts[genename], counts)]
                gene2ntx[genename] += 1

            # Column 2 of a collapsed row is the number of transcripts summed into it,
            # not a transcript ID; uncollapsed rows above keep their real ENST ID.
            for gene in order:
                counts_str = '\t'.join(str(round(m, 2)) for m in gene2counts[gene])
                outfh.write(f"{gene}\t{gene2ntx[gene]}\t{counts_str}\n")

        logger.info(f"Final counts ({level}-level, {len(order)} collapsed rows) "
                    f"written to {output_file}")

    # ===== Step 13: Reporting =====
    def _write_molecule_counts_custom_content(self, custom_dir):
        """Parse counts_molc_final.txt header lines into a MultiQC custom-content bargraph"""
        counts_file = os.path.join(self.config['output_dir'], 'counts_molc_final.txt')
        samples, unannot, annot = [], [], []
        with open(counts_file, 'r') as fh:
            for line in fh:
                if line.startswith('#samples'):
                    samples = line.strip().split('\t')[1:]
                elif line.startswith('#unannotatedmolc'):
                    unannot = line.strip().split('\t')[1:]
                elif line.startswith('#annotatedmolc'):
                    annot = line.strip().split('\t')[1:]
                    break

        molc_counts = {}
        outfile = os.path.join(custom_dir, 'smallrna_molecule_counts_mqc.tsv')
        with open(outfile, 'w') as fh:
            fh.write("# id: 'smallrna_molecule_counts'\n")
            fh.write("# section_name: 'SmallSeq Molecule Counts'\n")
            fh.write("# description: 'Annotated vs. unannotated small RNA molecule counts per sample after precursor removal and counting.'\n")
            fh.write("# plot_type: 'bargraph'\n")
            fh.write("# pconfig:\n")
            fh.write("#     id: 'smallrna_molecule_counts_plot'\n")
            fh.write("#     title: 'SmallSeq: Annotated vs Unannotated Molecule Counts'\n")
            fh.write("#     ylab: 'Molecule count'\n")
            fh.write("Sample\tannotated_molecules\tunannotated_molecules\n")
            for sample, a, u in zip(samples, annot, unannot):
                fh.write(f"{sample}\t{a}\t{u}\n")
                molc_counts[sample] = (a, u)

        return molc_counts

    def _write_umi_qc_custom_content(self, custom_dir):
        """Per-sample UMI health and saturation metrics.

        The Small-seq 5' adapter carries 8 'H' ribonucleotides (A/C/U only), so the UMI space
        is 3^8 = 6,561 — NOT 4^8. Two things follow, and neither is visible anywhere else:

          * A UMI containing G is impossible by design, so the observed G rate is a direct
            per-base error-rate readout for that sample.
          * With only 6,561 barcodes, an abundant small RNA can exhaust the UMI space at a
            single dedup group (locus x strand x read length). Molecule counts then saturate
            and undercount the most highly expressed species — worst in the deepest cells,
            which makes it a library-size-correlated artifact that normalisation won't remove.

        This reports the effect; it deliberately does not correct any counts.
        """
        UMI_SPACE = 3 ** 8  # 6,561

        rows = []
        outfile = os.path.join(custom_dir, 'smallrna_umi_qc_mqc.tsv')
        for sample in self.samples:
            fq = os.path.join(self.config['output_dir'], 'step1_umi_removed', sample,
                              f"{sample}_umiTrim.fq.gz")
            dedup_bam = os.path.join(self.config['output_dir'], 'step7_dedup', sample,
                                     f"{sample}_dedup.bam")
            if not os.path.exists(fq) or not os.path.exists(dedup_bam):
                logger.warning(f"Missing UMI QC inputs for {sample}, skipping in report")
                continue

            # --- UMI composition, from the tag umi_tools appended to each read name
            n_umi = n_with_g = g_bases = tot_bases = 0
            umi_counts = defaultdict(int)
            with gzip.open(fq, 'rt') as fh:
                for i, line in enumerate(fh):
                    if i % 4:
                        continue
                    umi = line.strip().rsplit('_', 1)[-1].split()[0]
                    if not umi:
                        continue
                    n_umi += 1
                    tot_bases += len(umi)
                    ng = umi.count('G')
                    g_bases += ng
                    if ng:
                        n_with_g += 1
                    # Only A/C/T UMIs are real. Anything with a G or an N is a miscall, and
                    # counting them would push distinct_UMIs above the 6,561 ceiling and
                    # inflate the diversity estimate.
                    if not ng and 'N' not in umi:
                        umi_counts[umi] += 1

            if not n_umi:
                continue

            # Effective diversity (inverse Simpson). Equals 6,561 only if every UMI is used
            # equally; T4 RNA ligase sequence bias makes real usage far more skewed, so
            # collisions -- and therefore undercounting -- happen sooner than nominal.
            tot = sum(umi_counts.values())
            u_eff = 1.0 / sum((c / tot) ** 2 for c in umi_counts.values()) if tot else 0.0

            # --- Saturation: molecules per dedup group. Post-dedup, 1 record == 1 molecule,
            # and each (locus, strand, length) group draws from its own independent UMI space.
            groups = defaultdict(int)
            bam = pysam.AlignmentFile(dedup_bam, "rb")
            for read in bam:
                if read.is_secondary or read.is_unmapped:
                    continue
                groups[(read.reference_id, read.reference_start,
                        read.is_reverse, read.query_length)] += 1
            bam.close()

            sizes = sorted(groups.values(), reverse=True) or [0]
            largest = sizes[0]
            over25 = sum(1 for s in sizes if s > 0.25 * UMI_SPACE)
            over50 = sum(1 for s in sizes if s > 0.50 * UMI_SPACE)
            saturated = "YES" if largest > 0.30 * UMI_SPACE else "no"
            if saturated == "YES":
                logger.warning(
                    f"{sample}: UMI saturation — largest dedup group holds {largest} molecules "
                    f"({100*largest/UMI_SPACE:.0f}% of the {UMI_SPACE}-UMI space). The most "
                    f"abundant species in this cell are undercounted.")

            rows.append((sample, 100*g_bases/tot_bases, 100*n_with_g/n_umi,
                         len(umi_counts), u_eff, largest,
                         100*largest/UMI_SPACE, over25, over50, saturated))

        with open(outfile, 'w') as fh:
            fh.write("# id: 'smallrna_umi_qc'\n")
            fh.write("# section_name: 'SmallSeq UMI QC'\n")
            fh.write("# description: 'UMI health and saturation. The 5-prime adapter uses 8 H bases "
                     "(A/C/U), so the UMI space is 3^8 = 6561 and a G in a UMI is always an error. "
                     "Dedup groups approaching that ceiling have undercounted molecules.'\n")
            fh.write("# plot_type: 'table'\n")
            fh.write("# pconfig:\n")
            fh.write("#     id: 'smallrna_umi_qc_table'\n")
            fh.write("#     title: 'SmallSeq: UMI QC and Saturation'\n")
            fh.write("Sample\tpct_G_bases\tpct_UMIs_with_G\tdistinct_UMIs\teff_UMI_diversity\t"
                     "largest_dedup_group\tpct_of_UMI_space\tgroups_over_25pct\tgroups_over_50pct\tsaturated\n")
            for r in rows:
                fh.write(f"{r[0]}\t{r[1]:.2f}\t{r[2]:.2f}\t{r[3]}\t{r[4]:.0f}\t"
                         f"{r[5]}\t{r[6]:.1f}\t{r[7]}\t{r[8]}\t{r[9]}\n")

        return {r[0]: r[1:] for r in rows}

    def _write_filtering_stats_custom_content(self, custom_dir):
        """Aggregate the step5/step8 per-sample stats files into a MultiQC custom-content bargraph"""
        def read_stats(stats_file):
            stats = {}
            with open(stats_file, 'r') as fh:
                for line in fh:
                    key, value = line.strip().split('\t')
                    stats[key] = value
            return stats

        filtering_stats = {}
        outfile = os.path.join(custom_dir, 'smallrna_filtering_stats_mqc.tsv')
        with open(outfile, 'w') as fh:
            fh.write("# id: 'smallrna_filtering_stats'\n")
            fh.write("# section_name: 'SmallSeq Filtering Stats'\n")
            fh.write("# description: 'Reads kept/removed by soft-clip removal (step 5) and precursor removal (step 8).'\n")
            fh.write("# plot_type: 'bargraph'\n")
            fh.write("# pconfig:\n")
            fh.write("#     id: 'smallrna_filtering_stats_plot'\n")
            fh.write("#     title: 'SmallSeq: Filtering Steps Read Counts'\n")
            fh.write("#     ylab: '# Reads'\n")
            fh.write("Sample\tstep5_kept\tstep5_removed\tstep8_kept\tstep8_removed\n")
            for sample in self.samples:
                softclip_file = os.path.join(self.config['output_dir'], 'step5_clipped_removed', sample, f"{sample}_softclip_stats.txt")
                precursor_file = os.path.join(self.config['output_dir'], 'step8_precursor_removed', sample, f"{sample}_precursor_stats.txt")
                if not os.path.exists(softclip_file) or not os.path.exists(precursor_file):
                    logger.warning(f"Missing filtering stats for {sample}, skipping in report")
                    continue

                softclip = read_stats(softclip_file)
                precursor = read_stats(precursor_file)
                step5_kept = int(softclip['kept_reads'])
                step5_removed = int(softclip['total_reads']) - step5_kept
                step8_kept = int(precursor['kept_reads'])
                step8_removed = int(precursor['total_reads']) - step8_kept
                fh.write(f"{sample}\t{step5_kept}\t{step5_removed}\t{step8_kept}\t{step8_removed}\n")
                filtering_stats[sample] = (step5_kept, step5_removed, step8_kept, step8_removed)

        return filtering_stats

    def step13_reporting(self):
        """Aggregate pipeline-specific stats into MultiQC custom content and generate the MultiQC report"""
        output_dir = self.config['output_dir']
        custom_dir = os.path.join(output_dir, 'multiqc_custom')
        self.safe_mkdir(custom_dir)

        molc_counts = self._write_molecule_counts_custom_content(custom_dir)
        filtering_stats = self._write_filtering_stats_custom_content(custom_dir)
        self._write_umi_qc_custom_content(custom_dir)

        cmd = f"multiqc {output_dir} -o {output_dir} -n multiqc_report"
        result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        if result.returncode != 0:
            logger.warning(f"MultiQC report generation warning: {result.stderr}")

        logger.info("Pipeline summary:")
        for sample in self.samples:
            if sample in molc_counts:
                annot, unannot = molc_counts[sample]
                logger.info(f"  {sample}: {annot} annotated / {unannot} unannotated molecules")
            if sample in filtering_stats:
                step5_kept, step5_removed, step8_kept, step8_removed = filtering_stats[sample]
                logger.info(f"  {sample}: soft-clip removal kept {step5_kept} (removed {step5_removed}), "
                            f"precursor removal kept {step8_kept} (removed {step8_removed})")
        logger.info(f"MultiQC report: {os.path.join(output_dir, 'multiqc_report.html')}")


def main():
    parser = argparse.ArgumentParser(
        description='Optimized SmallSeq Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Config file or individual parameters
    parser.add_argument('--config', help='YAML config file')
    parser.add_argument('--rawdata_dir', help='Raw data directory')
    parser.add_argument('--output_dir', help='Output directory')
    parser.add_argument('--genome_dir', help='Reference genome directory')
    parser.add_argument('--annotation', help='Gene annotation file (GenePred format)')
    parser.add_argument('--adapter_file', help='Adapter sequences file')
    parser.add_argument('--umi_pattern', default='NNNNNNNN', help='UMI pattern')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads')
    parser.add_argument('--max_read_len', type=int, default=40, help='Max read length for small RNA')
    parser.add_argument('--min_read_len', type=int, default=41, help='Min read length for precursor')
    parser.add_argument('--genome_fasta', help='Reference genome split fasta file directory')
    parser.add_argument('--legacy-count', action='store_true',
                        help='Use midpoint-based read assignment instead of whole-read overlap')
    parser.add_argument('--collapse-level', choices=['gene', 'transcript'], default='gene',
                        help='Report counts summed per gene (default), or per transcript with '
                             'only miRBase miRNAs collapsed (legacy behaviour)')
    parser.add_argument('--reset', action='store_true', help='Reset checkpoint and start from the beginning')
    parser.add_argument('--start-from', help='Start from a specific step (e.g., "STAR Alignment")')
    
    args = parser.parse_args()
    
    # Load config
    if args.config:
        if not YAML_AVAILABLE:
            raise ImportError("PyYAML required for config file support")
        with open(args.config, 'r') as f:
            config = yaml.safe_load(f)
    else:
        config = {
            'rawdata_dir': args.rawdata_dir,
            'genome_fasta': args.genome_fasta,
            'output_dir': args.output_dir,
            'genome_dir': args.genome_dir,
            'annotation': args.annotation,
            'threads': args.threads,
            'umi_pattern': args.umi_pattern,
            'max_read_len': args.max_read_len,
            'min_read_len': args.min_read_len,
            'legacy_count': args.legacy_count,
            'collapse_level': args.collapse_level,
        }
        if args.adapter_file:
            config['adapter_file'] = args.adapter_file
    
    # Run pipeline
    pipeline = SmallSeqPipeline(config)
    
    # Handle reset flag
    if args.reset:
        logger.info("Resetting checkpoint as requested")
        pipeline.reset_checkpoint()
    
    # Handle start-from flag
    if args.start_from:
        if args.start_from not in SmallSeqPipeline.STEP_NAMES:
            raise ValueError(
                f"Unknown step '{args.start_from}'. Valid steps: {SmallSeqPipeline.STEP_NAMES}"
            )
        idx = SmallSeqPipeline.STEP_NAMES.index(args.start_from)
        pipeline.completed_steps = set(SmallSeqPipeline.STEP_NAMES[:idx])
        pipeline.save_checkpoint()
        logger.info(f"Starting from step: {args.start_from}")
    
    pipeline.run()


if __name__ == '__main__':
    main()
