'''
Creation Date: 6/19/24
Updated by KyungTae Lee (postdoc, Ji Research Group)
Last update : 01/14/26

Description: Splits UMI and CB from 'read name' (string of text following #), discovers length of read, counts UMI-CB duplicates and collapses them. Next, it discovers how many reads have the same UMI, and sums the amount of reads with the same OMI. (Updated by KyungTae Lee) For dedupcliation, following criteria were used.

1. Remove reads whose mapping quality is 0 (unmapped)
2. Search for reads that share same Cell barcode - UMI
3. Group reads if they exhibit any overlap on genomic coordinates (on same chromosome)
4. Sort reads based on effective length(readlen X mapqual / max mapqual). Those reads with mapping quality equal to 255 are sorted by read length and placed in last priority. Select the first reads based on these criteria
5. If certain cell barcode-UMI come from three non-overlapping genomic loci, check for genomic loci that contains majority of reads (> 50%). If there is, select that loci and do the deduplication. Reads from any other loci are eliminated from downstream analysis.

Output: [OPTIONAL] <sample_name>.ref_span.txt - Noncollapsed version of reads
        [OPTIONAL] <sample_name>.collapsed.txt - Collapsed with counts
        <sample_name>.summary_output.tsv - BC, total reads with unique UMI, and total reads including duplicate UMIs.
        [OPTIONAL] BC_UMI.diff_genomic_loci.txt - cellBC_UMI, read_name, read_length, chrname, strand, ref_start, ref_end, map_qual of reads that contain same cell barcode-UMI and mapped to at least two different genomic loci (saved for later, newly added by KyungTae Lee)

Enivironment: Ran this script through long_reads - 'conda activate long_reads'
'''
import argparse
import os
import pysam
import csv
import sys
import random
import time
from collections import defaultdict

def parse_commandline():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam', '-b', help='input bam file', type=str, required=True)
    parser.add_argument("--statistics", action="store_true",
    help="If specified, output files for statistics of thee data")
    args = parser.parse_args()
#   print(args, file=sys.stderr)
    return args

def group_overlapping_reads(reads):
    """
    Groups reads by overlapping intervals -> sort by effective length and lenght if mapqual=255
    Args:
        reads (list of tuples): Each tuple contains (readname, read_length, read_start, read_end, chrname, strand, mapping_quality)
    
    Returns:
        list of lists: Each inner list is a group of overlapping reads.
    """
    # Sort reads by start position
    sorted_reads = sorted(reads, key=lambda x: x[2])  # Sort by read_start
    
    grouped_reads = []  # To store the groups
    current_group = []  # Current group of overlapping reads
    current_start, current_end = None, None  # Track current group's boundaries

    for read in sorted_reads:
        readname, read_length, read_start, read_end, chrname, strand, mapping_quality= read
    
        # If this is the first read or overlaps with the current group
        if not current_group or read_start <= current_end:
            current_group.append(read)
            # Update the group's end to the maximum end
            current_end = max(current_end, read_end) if current_end else read_end
        # If there is no overlap, sort reads and save them    
        else:
            # first put reads with mapq =255 at last and sort them based on read length
            read_255 = list(filter(lambda x: x[-1]==255, current_group))
            read_255= sorted(read_255, key=lambda x: -x[1])
            # second sort reads based on effective length(readlen X mapqual / max mapqual) 
            rest_read= list(filter(lambda x: x[-1] != 255, current_group))
            max_mapqual= max(list(map(lambda x: x[-1], rest_read)))
            rest_read= sorted(rest_read, key= lambda x: (float(x[1]) * float(x[-1])/ float(max_mapqual)), reverse=True)
            sorted_readgroup= rest_read+read_255
            # save the sorted read group and start a new one 
            grouped_reads.append(sorted_readgroup)
            current_group = [read]
            current_start, current_end = read_start, read_end

    # Add the last group
    if current_group:
        # first put reads with mapq =255 at last and sort them based on read length
        read_255 = list(filter(lambda x: x[-1]==255, current_group))
        read_255= sorted(read_255, key=lambda x: -x[1])
        # second sort reads based on effective length(readlen X mapqual / max mapqual) 
        rest_read= list(filter(lambda x: x[-1] != 255, current_group))
        max_mapqual= max(list(map(lambda x: x[-1], rest_read)))
        rest_read= sorted(rest_read, key= lambda x: (float(x[1]) * float(x[-1])/ float(max_mapqual)), reverse=True)
        sorted_readgroup= rest_read+read_255
        # save the last group
        grouped_reads.append(sorted_readgroup)

    return grouped_reads

def addBcCount(bc_countD, cell_barcode, dup_count, dedup_count=0):
    if not cell_barcode in bc_countD:
        bc_countD[cell_barcode]= {"dedup":0, "nodedup":0}
    # adding duplicated BC-UMI count
    bc_countD[cell_barcode]["nodedup"]+= dup_count 
    # adding de-duplicated BC-UMI count
    bc_countD[cell_barcode]["dedup"]+= dedup_count


def check_major_loci(outputlineL, outfile, bc_countD, cell_barcode):
    ## If there are groups of reads from at least 3 different genomic loci, check for major loci (50% > reads)
    total_readnum= sum(map(lambda x: float(x.rstrip().split("\t")[3]), outputlineL))
    outputlineL_sorted= sorted(outputlineL, key=lambda t: float(t.rstrip().split("\t")[3]), reverse= True)
    represent_line= outputlineL_sorted[0]
    max_readnum= float(outputlineL[0].rstrip().split("\t")[3])
    ratio= max_readnum/ float(total_readnum)
    if ratio > 0.5:
        addBcCount(bc_countD, cell_barcode, 0, dedup_count=1)
        outfile.write(represent_line)
    else:
        pass
   

def collapse_barcodes(barcode_dict, refspan_fn, diff_write):
    print("Processing Cell barcode (BC) - UMI collapsing")

    # Write updated data to a new file or overwrite original
    collapse_fn = refspan_fn.replace('.ref_span.txt', '.collapsed.txt')  
    outfile= open(collapse_fn, "w") 
    header= "\t".join(['Barcode_UMI', 'Representative Readname', 'Representative Length', 'Count'])+ "\n"
    outfile.write(header)
        
    # Collapsed read dictionary
    collapsed_readD= dict()
    # unique (dedup) and total (no dedup) read count per each cell barcodes dictionary
    bc_countD= dict()
    # barcode list containing multiple BC-UMI reads from different genomic loci
    bcumi_diffL= list()
    for bc_umi, readinfos in barcode_dict.items():
        cell_barcode= bc_umi.split("_")[0]
        # Remove unmapped read (mapping quality=0)
        readinfos= list(filter(lambda x: x[-1] != 0, readinfos))
    
        if len(readinfos)==0: # No mapped reads
            continue
        elif len(readinfos)==1: # Unique read with BC - UMI 
            readinfo= readinfos[0] ## since there is only read of this BC-UMI
            readname, read_length, reference_start, reference_end, chrname, strand, mapping_quality= readinfo
            # write unique BC-UMI read
            outputline= f"{bc_umi}\t{readname}\t{read_length}\t1\n"
            outfile.write(outputline)
            # add nodedup count information
            dup_count=1
            dedup_count=1
            addBcCount(bc_countD, cell_barcode, dup_count, dedup_count=dedup_count)        

        else: ## Multiple reads with same BC-UMI exist
            # Initialize dictionary to check genomic coordinate overlap
            collapse_readD= dict()
            # First, group reads based on chromosome name
            for i_readinfo in readinfos:
                readname, read_length, reference_start, reference_end, chrname, strand, mapping_quality= i_readinfo       
                if not chrname in collapse_readD:
                    collapse_readD[chrname]= list()
                collapse_readD[chrname].append(i_readinfo)
            # Initialize bcumi_diff count
            bcumi_diff=0
            outputlineL= list()
            # Second, group overlapping reads on same chromsome
            for i_chr in collapse_readD.keys():
                i_readinfos= collapse_readD[i_chr]
                if len(i_readinfos) ==1: # Unique BC-UMI in the same chromosome
                    # Add unique BC-UMI read info
                    readname, read_length, reference_start, reference_end, chrname, strand, mapping_quality= i_readinfos[0]
                    outputline= f"{bc_umi}\t{readname}\t{read_length}\t1\n"
                    outputlineL.append(outputline)
                    # Add dup count information
                    dup_count=1
                    addBcCount(bc_countD, cell_barcode, dup_count)        
                    bcumi_diff+=1
    
                else:  # Group reads if they show any overlap in genomic coordinates 
                    grouped_reads= group_overlapping_reads(i_readinfos)
                    for i_readgroup in grouped_reads:
                        represent_read= i_readgroup[0]
                        rep_readname= represent_read[0]
                        rep_readlen= represent_read[1]
                        read_count= len(i_readgroup)
                        # Add BC-UMI read info
                        outputline= f"{bc_umi}\t{rep_readname}\t{rep_readlen}\t{read_count}\n"
                        outputlineL.append(outputline)
                        # Add dup count information
                        dup_count=read_count
                        addBcCount(bc_countD, cell_barcode, dup_count)        
                        bcumi_diff+=1
            
            ## There are multiple BC-UMI combinations from different genomic loci
            if bcumi_diff >1: 
                bcumi_diffL.append(bc_umi)
            else:
                pass
            ## If there are same or more than 3 BC-UMI combinataions from different genomic loci, check if there is a loci with major reads (50% >=). 
            if bcumi_diff >= 3:
                check_major_loci(outputlineL, outfile, bc_countD, cell_barcode)
            else:
                for i_line in outputlineL:
                    outfile.write(i_line)
                    addBcCount(bc_countD, cell_barcode, 0, dedup_count=1)        
             
    outfile.close()
    print("")
    print(f"Collapsed barcodes and appended counts. Output saved to {collapse_fn}")

    if diff_write:
        bcumi_diff_output= "BC_UMI.diff_genomic_loci.txt"
        with open(bcumi_diff_output, "w") as outputopen:
            header= "cellBC_UMI\tread_name\tread_length\tchrname\tstrand\tref_start\tref_end\tmap_qual\n"
            outputopen.write(header)
            for bc_umi in bcumi_diffL:
                readinfos= barcode_dict[bc_umi]
                readinfos= list(filter(lambda x: x[-1] != 0, readinfos))
                for i_read in readinfos:
                    readname, read_length, reference_start, reference_end, chrname, strand, mapping_quality= i_read
                    outputline= f"{bc_umi}\t{readname}\t{read_length}\t{chrname}\t{strand}\t{reference_start}\t{reference_end}\t{mapping_quality}\n"
                    outputopen.write(outputline)
    else:
        pass
    return bc_countD, collapse_fn


def summarize(bc_countD, out_fn):
    print("Counting amount of reads in each UMI...")

    # Write results to output file
    output_file = out_fn.replace('ref_span.txt', 'summary_output.tsv')
    with open(output_file, 'w') as outfile:
        header= 'UMI\tUnique UMI reads for barcode \tTotal reads including dup UMI\n'
        outfile.write(header)
        for barcode in bc_countD.keys():
            unique_umi= bc_countD[barcode]["dedup"]
            total_umi= bc_countD[barcode]["nodedup"]
            outputline= f"{barcode}\t{unique_umi}\t{total_umi}\n"
            outfile.write(outputline)


    print(f"Summary written to {output_file}")


def main():
    current_time = time.localtime()
    formatted=time.strftime("%Y-%m-%d %H:%M:%S:", current_time)
    print(f"Current date and time: {formatted}")

    start_time = time.time()
    script_name = os.path.basename(__file__)
    print("Running", script_name)

    args = parse_commandline()

    bam_input = pysam.AlignmentFile(args.bam, 'rb')
    bam_basename = os.path.basename(args.bam)
    sample = bam_basename.split(".")[0]
    out_fn = sample + '.bc_umi.ref_span.txt'

    # Dictionary to store grouped entries by barcode
    barcode_dict = {}

    # Read bam file and group by unique format ofbarcode
    non_read=0
    for bamrd in bam_input.fetch(until_eof=True):
        query_parts = bamrd.query_name.split('#')
        barcode = query_parts[0]
        readparts = query_parts[1]
        readname = readparts
        # read_lenbgth= Number of aligned bases (M, =, X in CIGAR)
        try:
            read_length= sum(length for (operation, length) in bamrd.cigartuples if operation in (0, 7, 8))
        except TypeError:
            if bamrd.cigartuples is None and bamrd.is_unmapped:
                pass
            else:
                non_read+=1
            continue
        # Extract chrname
        chrname= bamrd.reference_name
        # Extract strand information
        strand = '-' if bamrd.is_reverse else '+'
        # Extract mapping quality
        map_qual= bamrd.mapping_quality

        # Groups reads by overarching barcode
        if barcode not in barcode_dict:
            barcode_dict[barcode] = []
        barcode_dict[barcode].append((readname, read_length, bamrd.reference_start, bamrd.reference_end, chrname, strand, map_qual))

    print (f"No Cigar reads that are not unmapped : {non_read}")
    # Write grouped and sorted entries to output file
    if args.statistics:
        with open(out_fn, 'w', newline='') as out_file:
            out_csv = csv.writer(out_file, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
            sorted_barcodes = sorted(barcode_dict.keys())
            for barcode in sorted_barcodes:
                barcode_dict[barcode].sort(key=lambda x: x[0])  # Sort each barcode/umi by readname
                for readname, ref_len, ref_start, ref_end, chrname, strand, map_qual in barcode_dict[barcode]:
                    out_csv.writerow([barcode, readname, ref_len, ref_start, ref_end, chrname, strand, map_qual])

    else:
        pass

    print("Processed and sorted, initializing barcode collapsing")

    # Call collapse function and summarize function
    bc_countD, collapse_fn= collapse_barcodes(barcode_dict, out_fn, args.statistics)

    summarize(bc_countD, out_fn)

    bam_input.close()

    end_time = time.time()
    elapsed=round(end_time - start_time, 3)
    print("Duration: ", elapsed, "seconds")

if __name__ == "__main__":
    main()
