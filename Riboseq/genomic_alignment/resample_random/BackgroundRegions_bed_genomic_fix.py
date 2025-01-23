##############################################################################
### BackgroundRegions_adapted_genomic.py takes a bedfile and a fasta as input and
### generates random regions of the fasta file with the length distribution as 
### given in the bed file and returns them in bed format
### for the genmic regions the selection of the random regions is made by the
### name, not by the chromsome (the chromsome was the name for transcriptomic
### alignment)
##############################################################################

#usage: python BackgroundRegions_adapted.py in.bed in.fasta out.bed
import sys
import random
from pybedtools import BedTool
import time


def get_length_dist(bed_file):
    #read in bedfile and obtain length distribution
    lengthdistribution=[]
    for line in bed_file:
        elems = line.split('\t')
        length=int(elems[2])-int(elems[1])
        lengthdistribution.append(length)
    return lengthdistribution


def get_random_bed_ids(lengthdistribution, UTR_bed_file):
    #select random keys for the fasta file, as many as there are 
        #entries in the bed file
    #do not allow for duplicate keys
        # Start time for Step 1
    three_prime_UTRs = BedTool(UTR_bed_file)
    random_regions_lengths = [(interval.name, interval.length) for interval in three_prime_UTRs]
    randomlist=[]
    for i in range(len(lengthdistribution)):
        required_length = lengthdistribution[i]
        filtered_intervals = [name for name, length in random_regions_lengths if length >= required_length]
        bed_names = list(set(filtered_intervals))
        random_key = random.choices(bed_names, k = 1)[0]
        #if key has already been sampled, randomly choose a new key
        assert random_key not in randomlist
        # while random_key in randomlist:
        #     random_key = random.choices(bed_names, k = 1)[0]
        randomlist.append(random_key)
        # filter the regions iteratively to prevent having the same name twice in the list
        random_regions_lengths = [(name, length) for name, length in random_regions_lengths if name != random_key]
    return randomlist, three_prime_UTRs

def write_bed_output(outname, randomlist, three_prime_UTRs, lengthdistribution):
    #for each length from the length distribution
        #get a random sample of the same length from the fasta sequence
    
    # random_list_reference = randomlist.copy()
    with open(outname, 'w') as out:
        i = len(randomlist)
        while i > 0:
            # start_time = time.time()
            length = lengthdistribution.pop()
            random_key = randomlist.pop()
            random_seq = three_prime_UTRs.filter(lambda interval: interval.name == random_key and interval.length >= length)[0]
            filter_bed_object_time = time.time()
            # print(f"Filter bed object for one region took {filter_bed_object_time - start_time:.2f} seconds.")
            # if len(random_seqs) == 0:
            #     raise ValueError(f"No intervals found for key: {random_key}")
            # elif len(random_seqs) == 1:
            #     random_seq = random_seqs[0]
            # else: 
            #     start_time = time.time()
            #     # if the region splits across exons, choose one of the intervals randomly
            #     # do not allow to sample several regions from the same UTR
            #     UTR_length = 0
            #     while length > UTR_length:
            #         interval_nr = random.randint(0, len(random_seqs) - 1)
            #         random_seq = random_seqs[interval_nr]
            #         UTR_length = random_seq.end - random_seq.start
            #     severl_regions_per_3prime_time = time.time()
            #     print(f"Several regions per 3 prime time took {severl_regions_per_3prime_time - start_time:.2f} seconds.")
            UTR_length = random_seq.end - random_seq.start
            assert length <= UTR_length
            if UTR_length == length: 
                start = random_seq.start
            else:
                start_time = time.time()
                start = random.randint(random_seq.start, UTR_length - length + random_seq.start)  
            end = start + length
            if random_seq.strand == '1':
                out.write(random_seq.chrom + '\t' + str(start) + '\t' + str(end) + '\t' + str(random_key)
                + '\t' + str(0) + '\t' + '+' + '\n')
            elif random_seq.strand == '-1':
                out.write(random_seq.chrom + '\t' + str(start) + '\t' + str(end) + '\t' + str(random_key)
                + '\t' + str(0) + '\t' + '-' + '\n')
            i -= 1 





if len(sys.argv) < 3:
        print('usage python BackgroundRegions_adapted_genomic.py in.bed UTR.bed out.bed')
else :
    random.seed(sys.argv[4])
    file = open(sys.argv[1],'r')
    lengthdistribution = get_length_dist(file)  
    
    randomlist, three_prime_UTRs = get_random_bed_ids(lengthdistribution, sys.argv[2])
    
    write_bed_output(sys.argv[3], randomlist, three_prime_UTRs, lengthdistribution)
    print(f'Random regions generated for iteration {sys.argv[4]}')
        
