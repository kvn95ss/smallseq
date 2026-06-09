#!/usr/bin/env python3
"""
Build Custom Annotation for SmallSeq Pipeline
Downloads and compiles standard transcript models (GENCODE), mature miRNAs (miRBase),
tRNAs (UCSC/GtRNAdb), and spike-ins into the extended GenePred format.

Usage:
    python3 build_annotation.py --output new_annots.gp --download
"""

import os
import sys
import argparse
import urllib.request
import gzip

# Default URLs for hg38/GRCh38 assembly
GENCODE_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.basic.annotation.gtf.gz"
MIRBASE_URL = "https://www.mirbase.org/download/hsa.gff3"
TRNA_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/tRNAs.txt.gz"

# Standard C. elegans spike-ins coordinates and structures as found in the original combined_annots.gp
SPIKE_INS = [
    ("P-cel-miR-41-5p", "P-cel-miR-41-5p", "+", 0, 24, 1),
    ("P-cel-miR-243-5p_2", "P-cel-miR-243-5p_2", "+", 0, 16, 2),
    ("P-cel-miR-60-5p", "P-cel-miR-60-5p", "+", 0, 21, 3),
    ("P-cel-miR-240-5p", "P-cel-miR-240-5p", "+", 0, 17, 4),
    ("P-cel-miR-59-3p", "P-cel-miR-59-3p", "+", 0, 20, 5),
    ("P-cel-miR-356a", "P-cel-miR-356a", "+", 0, 18, 6),
    ("P-cel-miR-792-5p", "P-cel-miR-792-5p", "+", 0, 25, 7),
    ("P-cel-miR-793", "P-cel-miR-793", "+", 0, 22, 8),
    ("P-cel-miR-1828", "P-cel-miR-1828", "+", 0, 19, 9),
    ("P-cel-miR-358-5p", "P-cel-miR-358-5p", "+", 0, 23, 10),
]

# Custom rRNA loci to append
CUSTOM_RRNA = [
    # Transcript ID, Gene Name, Chrom, Strand, Start, End
    ("rRNA45S", "rRNA45S", "chr21", "+", 8433222, 8446572),
    ("RNA5-8S5_ShortStack", "RNA5-8S5", "chr21", "+", 8395551, 8395789)
]

def download_file(url, dest):
    print(f"Downloading {url} to {dest}...")
    try:
        urllib.request.urlretrieve(url, dest)
        print("Download complete.")
    except Exception as e:
        print(f"Error downloading {url}: {e}")
        sys.exit(1)

def parse_attributes(attr_string):
    """Parse GTF attribute string into a dictionary"""
    attrs = {}
    for item in attr_string.strip().split(';'):
        if not item.strip():
            continue
        # Split by space (GTF standard)
        parts = item.strip().split(' ', 1)
        if len(parts) == 2:
            key, val = parts
            attrs[key] = val.replace('"', '').replace(';', '')
    return attrs

def parse_gencode_gtf(gtf_path):
    print(f"Parsing GENCODE GTF: {gtf_path}...")
    transcripts = []
    open_func = gzip.open if gtf_path.endswith('.gz') else open
    count = 0
    with open_func(gtf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            if parts[2] != 'transcript':
                continue
            
            chrom = parts[0]
            start = int(parts[3]) - 1  # Convert to 0-based
            end = int(parts[4])
            strand = parts[6]
            attrs = parse_attributes(parts[8])
            
            tx_id = attrs.get('transcript_id', '')
            gene_name = attrs.get('gene_name', tx_id)
            
            transcripts.append({
                'id': tx_id,
                'chrom': chrom,
                'strand': strand,
                'start': start,
                'end': end,
                'gene_name': gene_name
            })
            count += 1
    print(f"Loaded {count} transcripts from GENCODE.")
    return transcripts

def parse_mirbase_gff3(gff_path):
    print(f"Parsing miRBase GFF3: {gff_path}...")
    transcripts = []
    open_func = gzip.open if gff_path.endswith('.gz') else open
    count = 0
    with open_func(gff_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            # We want mature miRNAs (ID starts with MIMAT)
            if parts[2] not in ('miRNA', 'miRNA_primary_transcript'):
                continue
            
            attrs = {}
            for item in parts[8].strip().split(';'):
                if '=' in item:
                    k, v = item.split('=', 1)
                    attrs[k] = v
            
            tx_id = attrs.get('ID', '')
            if not tx_id.startswith('MIMAT'):  # Keep only mature miRNAs
                continue
            
            chrom = parts[0]
            start = int(parts[3]) - 1  # 0-based
            end = int(parts[4])
            strand = parts[6]
            gene_name = attrs.get('Name', tx_id)
            
            transcripts.append({
                'id': tx_id,
                'chrom': chrom,
                'strand': strand,
                'start': start,
                'end': end,
                'gene_name': gene_name
            })
            count += 1
    print(f"Loaded {count} mature miRNAs from miRBase.")
    return transcripts

def parse_ucsc_trnas(trna_path):
    print(f"Parsing UCSC tRNAs: {trna_path}...")
    transcripts = []
    open_func = gzip.open if trna_path.endswith('.gz') else open
    count = 0
    with open_func(trna_path, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 7:
                continue
            
            # UCSC database tRNAs track schema (has leading bin column):
            # 0: bin (e.g. 589)
            # 1: chrom (e.g. chr1)
            # 2: chromStart (e.g. 630995)
            # 3: chromEnd (e.g. 631061)
            # 4: name (e.g. nm-tRNA-Tyr-GTA-chr1-142)
            # 5: score (e.g. 1000)
            # 6: strand (e.g. -)
            chrom = parts[1]
            start = int(parts[2])
            end = int(parts[3])
            tx_id = parts[4]
            strand = parts[6]
            gene_name = parts[4]
            
            transcripts.append({
                'id': tx_id,
                'chrom': chrom,
                'strand': strand,
                'start': start,
                'end': end,
                'gene_name': gene_name
            })
            count += 1
    print(f"Loaded {count} tRNAs.")
    return transcripts

def write_extended_genepred(fh, tx_id, chrom, strand, start, end, gene_name, score=0):
    """Writes a single transcript model in extended GenePred format to file handle"""
    # Extended GenePred columns:
    # 0. bin (always 0 here)
    # 1. name (transcript_id)
    # 2. chrom
    # 3. strand
    # 4. txStart
    # 5. txEnd
    # 6. cdsStart (same as txStart)
    # 7. cdsEnd (same as txEnd)
    # 8. exonCount (1 exon representing the full locus span)
    # 9. exonStarts (start,)
    # 10. exonEnds (end,)
    # 11. score (default 0)
    # 12. name2 (gene name/symbol)
    # 13. cdsStartStat (none)
    # 14. cdsEndStat (none)
    # 15. exonFrames (-1)
    line_parts = [
        "0",
        tx_id,
        chrom,
        strand,
        str(start),
        str(end),
        str(start),
        str(end),
        "1",
        f"{start},",
        f"{end},",
        str(score),
        gene_name,
        "none",
        "none",
        "-1"
    ]
    fh.write("\t".join(line_parts) + "\n")

def main():
    parser = argparse.ArgumentParser(description="Generate Custom extended GenePred Annotation for SmallSeq Pipeline")
    parser.add_argument("-g", "--gencode", help="Path to GENCODE GTF file")
    parser.add_argument("-m", "--mirbase", help="Path to miRBase GFF3 file")
    parser.add_argument("-t", "--trna", help="Path to UCSC tRNAs txt.gz file")
    parser.add_argument("-o", "--output", default="combined_annots.gp", help="Output GenePred file path")
    parser.add_argument("--download", action="store_true", help="Download missing files automatically")
    args = parser.parse_args()

    # Create directories if they don't exist
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    gencode_file = args.gencode
    mirbase_file = args.mirbase
    trna_file = args.trna

    # Handle automatic downloads
    if args.download:
        if not gencode_file:
            gencode_file = "gencode.v46.basic.annotation.gtf.gz"
            if not os.path.exists(gencode_file):
                download_file(GENCODE_URL, gencode_file)
        if not mirbase_file:
            mirbase_file = "hsa.gff3"
            if not os.path.exists(mirbase_file):
                download_file(MIRBASE_URL, mirbase_file)
        if not trna_file:
            trna_file = "tRNAs.txt.gz"
            if not os.path.exists(trna_file):
                download_file(TRNA_URL, trna_file)

    # Validate that we have inputs
    if not gencode_file or not os.path.exists(gencode_file):
        print("Error: GENCODE GTF file not found. Use --gencode or run with --download.")
        sys.exit(1)
    if not mirbase_file or not os.path.exists(mirbase_file):
        print("Error: miRBase GFF3 file not found. Use --mirbase or run with --download.")
        sys.exit(1)
    if not trna_file or not os.path.exists(trna_file):
        print("Error: tRNA file not found. Use --trna or run with --download.")
        sys.exit(1)

    # Parse all databases
    gencode_transcripts = parse_gencode_gtf(gencode_file)
    mirbase_miRNAs = parse_mirbase_gff3(mirbase_file)
    trnas = parse_ucsc_trnas(trna_file)

    # Write combined annotation file
    print(f"Writing combined GenePred file: {args.output}...")
    with open(args.output, "w") as out:
        # 1. Write miRNAs
        for tx in mirbase_miRNAs:
            write_extended_genepred(out, tx['id'], tx['chrom'], tx['strand'], tx['start'], tx['end'], tx['gene_name'])
        
        # 2. Write GENCODE transcripts
        for tx in gencode_transcripts:
            write_extended_genepred(out, tx['id'], tx['chrom'], tx['strand'], tx['start'], tx['end'], tx['gene_name'])
        
        # 3. Write tRNAs
        for tx in trnas:
            write_extended_genepred(out, tx['id'], tx['chrom'], tx['strand'], tx['start'], tx['end'], tx['gene_name'])

        # 4. Write Custom rRNA
        for tx_id, gene_name, chrom, strand, start, end in CUSTOM_RRNA:
            write_extended_genepred(out, tx_id, chrom, strand, start, end, gene_name)

        # 5. Write Spike-ins
        for tx_id, gene_name, strand, start, end, score in SPIKE_INS:
            # Note: Spike-ins map to their own artificial chromosome named after the spike-in
            write_extended_genepred(out, tx_id, tx_id, strand, start, end, gene_name, score=score)

    print("Success! Custom annotation file generated successfully.")

if __name__ == "__main__":
    main()
