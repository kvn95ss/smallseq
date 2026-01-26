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


class SmallSeqPipeline:
    """Main pipeline class for SmallSeq data processing"""
    
    def __init__(self, config):
        self.config = config
        self.samples = []
        self.validate_config()
        
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
        self.config.setdefault('max_read_len', 40)
        self.config.setdefault('min_read_len', 41)
        self.config.setdefault('adapter_file', '../adapters/cutadapt_3prime.fa')
        self.config.setdefault('allowed_5p_clip', 0)
        self.config.setdefault('allowed_3p_clip', 3)
        self.config.setdefault('dedup_method', 'adjacency')
        
    def safe_mkdir(self, path):
        """Create directory if it doesn't exist"""
        if not os.path.exists(path):
            os.makedirs(path, mode=0o774)
            logger.info(f"Created directory: {path}")
    
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
        steps = [
            ("UMI Removal", self.step1_remove_umi),
            ("Adapter Trimming", self.step2_trim_adapters),
            ("STAR Alignment", self.step3_star_alignment),
#            ("SAM Processing", self.step4_process_sam),
            ("Soft-clip Removal", self.step5_remove_softclipped),
            ("Read Length Filtering", self.step6_filter_by_length),
            ("UMI Deduplication", self.step7_umi_dedup),
            ("Precursor Removal", self.step8_remove_precursors),
            ("Count Generation", self.step9_count_smallrnas),
            ("Count Merging", self.step10_merge_counts),
            ("miRNA Collapsing", self.step11_collapse_mirnas)
        ]
        
        for step_name, step_func in steps:
            logger.info(f"\n{'='*60}")
            logger.info(f"Step: {step_name}")
            logger.info(f"{'='*60}")
            try:
                step_func()
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
        
        if not os.path.exists(input_fq):
            logger.warning(f"Input file not found: {input_fq}")
            return
        
        cmd = f"cutadapt -a file:{self.config['adapter_file']} " \
              f"-e 0.1 -O 1 -u 2 --quiet --minimum-length 18 " \
              f"-o {output_fq} {input_fq}"
        
        #subprocess.run(cmd, shell=True, check=True)
        result = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
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
        
        fq_list = os.path.join(output_dir, "star_fastqs.txt")
        rg_list = os.path.join(output_dir, "star_readgroups.txt")
        samples_used = []
        
        # 1) Collect FASTQs + read groups
        with open(fq_list, "w") as fq_fh, open(rg_list, "w") as rg_fh:
            for sample in self.samples:
                fq = os.path.join(input_dir, sample, f"{sample}.fastq.gz")
                if not os.path.exists(fq):
                    logger.warning(f"Input file not found: {fq}")
                    continue
                fq_fh.write(f"{fq}\t-\tID:{sample}\tSM:{sample}\tPL:ILLUMINA\n")
                rg_fh.write(f"ID:{sample} SM:{sample} PL:ILLUMINA\n")
                samples_used.append(sample)
        
        if not samples_used:
            raise RuntimeError("No valid FASTQs found for STAR alignment")
        
        # 2) Run STAR
        prefix = os.path.join(output_dir, "all_samples_")
        star_cmd = (
            f"STAR "
            f"--runThreadN {self.config['threads']} "
            f"--genomeDir {self.config['genome_dir']} "
            f"--readFilesManifest {fq_list} "
            f"--readFilesCommand zcat "
            f"--outSAMtype BAM SortedByCoordinate "
            f"--outSAMstrandField intronMotif "
            f"--outFilterMultimapNmax 50 "
            f"--outFilterScoreMinOverLread 0 "
            f"--outFilterMatchNmin 18 "
            f"--outFilterMatchNminOverLread 0 "
            f"--outFilterMismatchNoverLmax 0.04 "
            f"--alignIntronMax 1 "
            f"--outFileNamePrefix {prefix} "
            f"--outSAMattributes NH HI AS nM RG"
        )
        
        result = subprocess.run(
            star_cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        if result.returncode != 0:
            raise RuntimeError(f"STAR failed:\n{result.stderr}")
        
        # 3) Split BAM by read group
        combined_bam = os.path.join(
            output_dir, "all_samples_Aligned.sortedByCoord.out.bam"
        )
        if not os.path.exists(combined_bam):
            raise FileNotFoundError("Combined STAR BAM not found")
        
        split_cmd = f"samtools split -@ {self.config['threads']} -M -1 -d RG {combined_bam}"
        result = subprocess.run(
            split_cmd,
            shell=True,
            cwd=output_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        if result.returncode != 0:
            raise RuntimeError(f"samtools split failed:\n{result.stderr}")
        
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
        
        inbam = os.path.join(input_dir, sample, f"{sample}.bam")
        outbam = os.path.join(sample_out, f"{sample}.bam")
        
        cmd = f"samtools view -h {inbam} | " \
              f"awk 'length($10) <= {self.config['max_read_len']} || $1~\"@\"' | " \
              f"samtools view -bS - > {outbam}"
        
        #subprocess.run(cmd, shell=True, check=True)
        result = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
        if result.returncode !=0:
            logger.warning(f"Filtering warning for {sample}: {result.stderr}")
        
        pysam.index(outbam)
    
    def step6_filter_by_length(self):
        """Filter reads by length"""
        output_dir = os.path.join(self.config['output_dir'], f'step6_max{self.config["max_read_len"]}')
        self.safe_mkdir(output_dir)
        
        with Pool(self.config['threads']) as pool:
            pool.map(self._process_readlen_sample, self.samples)
    
    # ===== Step 7: UMI Deduplication =====
    def _process_dedup_sample(self, sample):
        """Process a single sample for UMI deduplication"""
        input_dir = os.path.join(self.config['output_dir'], f'step6_max{self.config["max_read_len"]}')
        output_dir = os.path.join(self.config['output_dir'], 'step7_dedup')
        
        sample_out = os.path.join(output_dir, sample)
        self.safe_mkdir(sample_out)
        
        inbam = os.path.join(input_dir, sample, f"{sample}.bam")
        outbam = os.path.join(sample_out, f"{sample}.bam")
        logfile = os.path.join(sample_out, "dedup.log")
        
        cmd = f"umi_tools dedup --method {self.config['dedup_method']} " \
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
    def _process_precursor_sample(self, sample, gf):
        """Process a single sample for precursor removal"""
        input_dir = os.path.join(self.config['output_dir'], 'step7_dedup')
        output_dir = os.path.join(self.config['output_dir'], 'step8_precursor_removed')
        
        sample_out = os.path.join(output_dir, sample)
        self.safe_mkdir(sample_out)
        
        inbam = os.path.join(input_dir, sample, f"{sample}.bam")
        outbam_tmp = os.path.join(sample_out, f"{sample}_tmp.bam")
        outbam = os.path.join(sample_out, f"{sample}.bam")
        
        inbam_obj = pysam.AlignmentFile(inbam, "rb")
        outbam_obj = pysam.AlignmentFile(outbam_tmp, "wb", template=inbam_obj)
        
        for read in inbam_obj:
            readlen = len(read.seq)
            
            if readlen <= 35:  # Keep short reads
                outbam_obj.write(read)
                continue
            
            readchr = inbam_obj.get_reference_name(read.reference_id)
            readstart = read.pos + 1
            readend = read.reference_end
            
            minRlen = self.config['min_read_len']
            upperlimit = minRlen - readlen
            
            keep_read = True
            if read.flag == 0:  # Forward strand
                bpwindow = gf.get_seq_from_to(readchr, readend+1, readend+upperlimit)
                patterns = {
                    minRlen-1: ["T", "A"],
                    minRlen-2: ["TG", "AA"],
                    minRlen-3: ["TGG", "AAA"],
                    minRlen-4: ["TGGA", "AAAA"],
                    minRlen-5: ["TGGAA", "AAAAA"]
                }
                if readlen in patterns and bpwindow in patterns[readlen]:
                    keep_read = False
                    
            elif read.flag == 16:  # Reverse strand
                bpwindow = gf.get_seq_from_to(readchr, readstart-upperlimit, readstart-1)
                patterns = {
                    minRlen-1: ["A", "T"],
                    minRlen-2: ["CA", "TT"],
                    minRlen-3: ["CCA", "TTT"],
                    minRlen-4: ["TCCA", "TTTT"],
                    minRlen-5: ["TTCCA", "TTTTT"]
                }
                if readlen in patterns and bpwindow in patterns[readlen]:
                    keep_read = False
            
            if keep_read:
                outbam_obj.write(read)
        
        outbam_obj.close()
        inbam_obj.close()
        
        # Sort and index
        pysam.sort('-o', outbam, outbam_tmp)
        pysam.index(outbam)
        os.remove(outbam_tmp)
    
    def step8_remove_precursors(self):
        """Remove reads from precursor RNAs based on genomic context"""
        output_dir = os.path.join(self.config['output_dir'], 'step8_precursor_removed')
        self.safe_mkdir(output_dir)
        
        # Import GenomeFetch
        sys.path.insert(0, os.path.dirname(__file__))
        try:
            from GenomeFetch import GenomeFetch
        except ImportError:
            logger.error("GenomeFetch.py not found. Please ensure it's in the same directory.")
            raise
        
        gf = GenomeFetch(genomedir=self.config['genome_fasta'])
        
        # Process samples sequentially due to GenomeFetch not being thread-safe
        for sample in self.samples:
            self._process_precursor_sample(sample, gf)
    
    # ===== Step 9: Count Small RNAs =====
    def _process_count_sample(self, sample, intervals, coord2geneid, geneid2name, geneidlist):
        """Process a single sample for counting"""
        input_dir = os.path.join(self.config['output_dir'], 'step8_precursor_removed')
        output_dir = os.path.join(self.config['output_dir'], 'step9_counts')
        
        sample_out = os.path.join(output_dir, sample)
        self.safe_mkdir(sample_out)
        
        inbam = os.path.join(input_dir, sample, f"{sample}.bam")
        outfile = os.path.join(sample_out, f"{sample}_Count.txt")
        
        def find_overlaps(chrom, pos, strand):
            """Find genes overlapping a position"""
            overlaps = []
            if chrom in intervals and strand in intervals[chrom]:
                for start, end, geneid in intervals[chrom][strand]:
                    if start <= pos <= end:
                        overlaps.append(f"{chrom}:{start+1}-{end}:{strand}")
            return overlaps
        
        inbam_obj = pysam.AlignmentFile(inbam, "rb")
        
        read2overlaps = defaultdict(list)
        
        for read in inbam_obj:
            readchr = inbam_obj.get_reference_name(read.reference_id)
            readstart = read.pos
            readend = read.reference_end
            strand = "-" if read.is_reverse else "+"
            
            midpos = (readstart + readend) // 2
            overlaps = find_overlaps(readchr, midpos, strand)
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
                        for coord in overlaps:
                            geneid = coord2geneid.get(coord, 'NA')
                            if geneid not in geneid2counts:
                                geneid2counts[geneid] = 0
                            geneid2counts[geneid] += 1 / annot_count
            else:
                num_unannot += 1
        
        num_annot = sum(v for k, v in geneid2counts.items() if k != 'NA' and not k.startswith('P-cel'))
        
        # Write output
        with open(outfile, 'w') as fh:
            fh.write(f"#samples\t{sample}\n")
            fh.write(f"#unannotatedmolc\t{num_unannot}\n")
            fh.write(f"#annotatedmolc\t{num_annot}\n")
            for geneid in geneidlist:
                fh.write(f"{geneid2name[geneid]}\t{geneid}\t{geneid2counts.get(geneid, 0)}\n")
    
    def step9_count_smallrnas(self):
        """Count small RNAs using annotation"""
        output_dir = os.path.join(self.config['output_dir'], 'step9_counts')
        self.safe_mkdir(output_dir)
        
        # Load annotation
        geneid2name = {}
        coord2geneid = {}
        geneidlist = []
        
        # Simple interval tree implementation
        intervals = defaultdict(lambda: defaultdict(list))
        
        for line in open(self.config['annotation'], 'r'):
            p = line.split()
            chrom, start, end, strand, geneid, genename = p[2], int(p[4]), int(p[5]), p[3], p[1], p[12]
            
            intervals[chrom][strand].append((start, end, geneid))
            coord = f"{chrom}:{start+1}-{end}:{strand}"
            coord2geneid[coord] = geneid
            geneid2name[geneid] = genename
            geneidlist.append(geneid)
        
        # Process samples sequentially to avoid shared state issues
        for sample in self.samples:
            self._process_count_sample(sample, intervals, coord2geneid, geneid2name, geneidlist)
    
    # ===== Step 10: Merge Counts =====
    def step10_merge_counts(self):
        """Merge count files from all samples"""
        input_dir = os.path.join(self.config['output_dir'], 'step9_counts')
        output_file = os.path.join(self.config['output_dir'], 'counts_molc.txt')
        
        count_files = []
        for sample in self.samples:
            count_file = os.path.join(input_dir, sample, f"{sample}_Count.txt")
            if os.path.exists(count_file):
                count_files.append(count_file)
        
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
    
    # ===== Step 11: Collapse miRNAs =====
    def step11_collapse_mirnas(self):
        """Collapse multi-mapping miRNAs"""
        input_file = os.path.join(self.config['output_dir'], 'counts_molc.txt')
        output_file = os.path.join(self.config['output_dir'], 'counts_molc_final.txt')
        
        gene2molc = {}
        gene2trid = {}
        
        with open(output_file, 'w') as outfh:
            for line in open(input_file, 'r'):
                if line.startswith('#'):
                    outfh.write(line)
                else:
                    p = line.strip().split('\t')
                    genename = p[0]
                    trans_ids = p[1]
                    gene2trid[genename] = trans_ids
                    
                    if genename.startswith("hsa"):  # miRBase miRNAs
                        molc_counts = list(map(float, p[2:]))
                        zeros = [0] * len(molc_counts)
                        gene2molc[genename] = [i+j for i, j in zip(gene2molc.get(genename, zeros), molc_counts)]
                    else:
                        outfh.write(line)
            
            # Write collapsed miRNAs
            for gene in gene2molc:
                counts_str = '\t'.join([str(round(m, 2)) for m in gene2molc[gene]])
                outfh.write(f"{gene}\t{gene2trid[gene]}\t{counts_str}\n")
        
        logger.info(f"Final counts written to {output_file}")
        
        # Clean up intermediate file
        os.remove(input_file)


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
        }
        if args.adapter_file:
            config['adapter_file'] = args.adapter_file
    
    # Run pipeline
    pipeline = SmallSeqPipeline(config)
    pipeline.run()


if __name__ == '__main__':
    main()
