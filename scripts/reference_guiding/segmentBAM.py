#!/usr/bin/env python

import sys
import os
import math
import pysam
import argparse

import multiprocessing as mp

from statistics import mean

def printBinInfo(info):

	if (info[3] + info[4] != 0): print(info[0] + '\t' + info[1] + '\t' +  info[2] + '\t' + str(info[3]) + '\t' + str(info[4]))

def outputUnmappedLowQualLR(out_prefix_filename, mapq, ref_unbin, multithreaded):

	global lr_bams

	nb_base_lr = 0

	out_f = open(out_prefix_filename + "_lr_unknown.fq", "w") # Output file

	for bamf in lr_bams:

		it_bamf = bamf.fetch(multiple_iterators=multithreaded, until_eof=True) # fetch() also gets unmapped reads here

		for record in it_bamf:

			# Discard 0-length reads, unmapped reads, secondary and supplementary alignments.
			if (record.query_length != 0) and ((record.is_unmapped == True) or ((record.is_secondary == False) and (record.is_supplementary == False) and ((record.mapping_quality < mapq) or (record.reference_name in ref_unbin)))):

				# Output LR
				out_f.write('>' + record.query_name + '\n' + record.query_sequence + '\n')

				nb_base_lr += record.query_length
				
				# If has qualities, just output those
				if (record.query_qualities != None):
				
					qual = [chr(x+33) for x in record.query_qualities]
				
					out_f.write('+' + '\n' + "".join(qual) + '\n')

	out_f.close()

	return ("unknown_lr", os.path.abspath(out_prefix_filename + "_lr_unknown.fq"), "NA", nb_base_lr, 0)

def outputUnmappedLowQualSR(out_prefix_filename, mapq, qs_read, qs_base, multithreaded):

	global sr_bams

	nb_base_sr = 0

	out_f = open(out_prefix_filename + "_sr_unmapped.fa", "w") # Output file for unmapped or low qual SR
	out_junk_f = open(out_prefix_filename + "_sr_discarded.fa", "w") # Ouput file for discarded unmapped or low qual SR (junk)

	for bamf in sr_bams:

		it_bamf = bamf.fetch(multiple_iterators=multithreaded, until_eof=True) # fetch() also gets unmapped reads here

		for record in it_bamf:

			# Discard 0-length reads, unmapped reads, secondary and supplementary alignments.
			if (record.query_length != 0) and ((record.is_unmapped == True) or ((record.is_secondary == False) and (record.is_supplementary == False) and (record.mapping_quality < mapq))):

				nb_base_sr += record.query_length

				if (record.query_qualities == None): out_f.write('>' + record.query_name + '\n' + record.query_sequence + '\n')
				elif (mean(record.query_qualities) < qs_read):

					out_junk_f.write('>' + record.query_name + '\n' + record.query_sequence + '\n')

					nb_base_sr -= record.query_length

				elif (qs_base != 0):

					pos_low_qual = [ i for i in range(len(record.query_qualities)) if (record.query_qualities[i] < qs_base) ]
					query_sequence = list(record.query_sequence)
					
					for pos in pos_low_qual: query_sequence[pos] = 'N'

					out_f.write('>' + record.query_name + '\n' + ''.join(query_sequence) + '\n')

				else: out_f.write('>' + record.query_name + '\n' + record.query_sequence + '\n')



	out_f.close()
	out_junk_f.close()

	return ("unmapped_sr", os.path.abspath(out_prefix_filename + "_sr_unmapped.fa"), "NA", nb_base_sr, 0)
	
def segmentBAM(out_prefix_filename, chr_name, lr_start_pos_ref, lr_end_pos_ref, chr_len, mapq_sr, mapq_lr, qs_sr, multithreaded):

	global sr_bams
	global lr_bams

	buffer_sz = 1000000 # Size buffer to add to SR reads boundaries

	sr_start_pos_ref = chr_len # SR start position
	sr_end_pos_ref = 0 # SR end position

	bin_name = chr_name + "_" + str(lr_start_pos_ref) # Bin name

	out_sr_filename = out_prefix_filename + "_sr_" + bin_name + ".fa" # Output name SR file
	out_lr_filename = out_prefix_filename + "_lr_" + bin_name + ".fq" # Output name LR file

	out_sr_f = open(out_sr_filename, "w") # Output file for SR
	out_lr_f = open(out_lr_filename, "w") # Output file for LR

	nb_base_sr = 0 # Number of SR extracted for that region
	nb_base_lr = 0 # Number of LR extracted for that region

	for lr_bamf in lr_bams:

		it_lr_bamf = lr_bamf.fetch(chr_name, lr_start_pos_ref, lr_end_pos_ref, multiple_iterators=multithreaded)

		for record in it_lr_bamf:

			# If long read starts in the region and its mapq is good enough
			if (record.reference_start >= lr_start_pos_ref) and (record.reference_start <= lr_end_pos_ref) and (record.mapping_quality >= mapq_lr):

				# Discard 0-length reads, unmapped reads, secondary and supplementary alignments.
				if (record.query_length != 0) and (record.is_unmapped == False) and (record.is_secondary == False) and (record.is_supplementary == False):

					# Set the boundaries where we need to fetch the SR
					sr_start_pos_ref = min(sr_start_pos_ref, record.reference_start)
					sr_end_pos_ref = max(sr_end_pos_ref, record.reference_start + record.reference_length - 1)

					# Output LR
					out_lr_f.write('>' + record.query_name + '\n' + record.query_sequence + '\n')

					# If has qualities, just output those
					if (record.query_qualities != None):

						qual = [chr(x+33) for x in record.query_qualities]

						out_lr_f.write('+' + '\n' + "".join(qual) + '\n')

					nb_base_lr += record.query_length # Increase count of LR extracted

	if (nb_base_lr != 0): # If there are long reads to correct for that region

		# Increase SR boundaries by buffer_sz
		sr_start_pos_ref = max(sr_start_pos_ref - buffer_sz, 0)
		sr_end_pos_ref = min(sr_end_pos_ref + buffer_sz, chr_len)

		for sr_bamf in sr_bams:

			it_sr_bamf = sr_bamf.fetch(chr_name, sr_start_pos_ref, sr_end_pos_ref, multiple_iterators=multithreaded)

			for record in it_sr_bamf:

				# If short read starts in the region and its mapq is good enough
				if (record.reference_start >= sr_start_pos_ref) and (record.reference_start <= sr_end_pos_ref) and (record.mapping_quality >= mapq_sr):

					# Discard 0-length reads, unmapped reads, secondary and supplementary alignments.
					if (record.query_length != 0) and (record.is_unmapped == False) and (record.is_secondary == False) and (record.is_supplementary == False):

						if (qs_sr != 0) and (record.query_qualities != None):

							pos_low_qual = [ i for i in range(len(record.query_qualities)) if (record.query_qualities[i] < qs_sr) ]
							query_sequence = list(record.query_sequence)
							
							for pos in pos_low_qual: query_sequence[pos] = 'N'

							out_sr_f.write('>' + record.query_name + '\n' + ''.join(query_sequence) + '\n')

						# Output LR
						else: out_sr_f.write('>' + record.query_name + '\n' + record.query_sequence + '\n')

						nb_base_sr += record.query_length # Increase count of LR extracted

	out_sr_f.close()
	out_lr_f.close()

	if (nb_base_lr == 0) and (nb_base_sr == 0):

		os.remove(os.path.abspath(out_sr_filename))
		os.remove(os.path.abspath(out_lr_filename))

		return (bin_name, "NA", "NA", 0, 0)

	else: return (bin_name, os.path.abspath(out_sr_filename), os.path.abspath(out_lr_filename), nb_base_sr, nb_base_lr)

def checkReferenceCompatibility(sr_filenames, lr_filenames, force_inter_ref):

	ref_s_inter = set()
	ref_s_union = set()

	ref_d_inter = {}
	ref_d_diff = {}

	filenames = sr_filenames + lr_filenames
	first_file = True

	for filename in filenames:

		ref_file = set()
		bamf = pysam.AlignmentFile(filename, "rb")

		for chr_name in bamf.references: ref_file.add((chr_name, bamf.get_reference_length(chr_name)))

		if (first_file == True):

			ref_s_inter = ref_file
			ref_s_union = ref_file
			first_file = False

		else:

			ref_s_inter = ref_s_inter.intersection(ref_file)
			ref_s_union = ref_s_union.union(ref_file)

	ref_s_diff = ref_s_union.difference(ref_s_inter)

	if (len(ref_s_diff) != 0) and (force_inter_ref == False):

		sys.exit("Input BAM files have different reference chromsomes/contigs. Use option --intersection_ref to force using the intersection.")

	for contig in ref_s_inter: ref_d_inter[contig[0]] = contig[1]
	for contig in ref_s_diff: ref_d_diff[contig[0]] = contig[1]

	return (ref_d_inter, ref_s_diff)

if __name__ == '__main__':

	# Default values
	sr_filenames = []
	lr_filename = []
	out_prefix_filename = ""
	nb_threads = 1
	len_segment = 5000000
	mapq_sr = 0
	mapq_lr = 30
	qs_sr = 10
	qs_unmap_sr = 20
	force_inter_ref = False

	# Parse arguments
	parser = argparse.ArgumentParser(prog='segmentBAM', description='Segment BAM files of short reads and long reads', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	required = parser.add_argument_group('Required arguments')

	required.add_argument('-s', '--short_read_bam', action='append', help='Filename of a short read bam file', required=True)
	required.add_argument('-l', '--long_read_bam', action='append', help='Filename of a long read bam file', required=True)
	required.add_argument('-o', '--out_prefix_filename', action='store', help='Prefix of the output filenames', required=True)

	optional = parser.add_argument_group('Optional arguments')

	optional.add_argument('-t', '--threads', action='store', help='Number of threads to use', default=nb_threads, required=False)
	optional.add_argument('-b', '--buffer_size', action='store', help='Length of segments in bp', default=len_segment, required=False)
	optional.add_argument('-m', '--mapq_short', action='store', help='Minimum MAPQ of short reads', default=mapq_sr, required=False)
	optional.add_argument('-n', '--mapq_long', action='store', help='Minimum MAPQ of long reads', default=mapq_lr, required=False)
	optional.add_argument('-q', '--qs_sr', action='store', help='Minimum quality score of short read bases', default=qs_sr, required=False)
	optional.add_argument('-u', '--qs_unmap_sr', action='store', help='Minimum mean quality score of unmapped short reads', default=qs_unmap_sr, required=False)

	optional.add_argument('--intersection_ref', action='store_true', help='Force using the intersection of reference contigs if input BAMs have different references.', required=False)

	args = parser.parse_args()
	args_d = vars(args)

	for k,v in args_d.items():

		if ((k == "s") or (k == "short_read_bam")): sr_filenames = v
		elif ((k == "l") or (k == "long_read_bam")): lr_filenames = v
		elif ((k == "t") or (k == "threads")): nb_threads = int(v)
		elif ((k == "o") or (k == "out_prefix_filename")): out_prefix_filename = v
		elif ((k == "b") or (k == "buffer_size")): len_segment = int(v)
		elif ((k == "m") or (k == "mapq_short")): mapq_sr = int(v)
		elif ((k == "n") or (k == "mapq_long")): mapq_lr = int(v)
		elif ((k == "q") or (k == "qs_sr")): qs_sr = int(v)
		elif ((k == "u") or (k == "qs_unmap_sr")): qs_unmap_sr = int(v)
		elif (k == "intersection_ref"): force_inter_ref = True

	# Minor consistency check on arguments
	if (len(sr_filenames) == 0): sys.exit("No input short read BAM provided as input")
	if (len(lr_filenames) == 0): sys.exit("No input long read BAM provided as input")
	if (nb_threads <= 0): sys.exit("Cannot use less than 1 thread")
	if (nb_threads > mp.cpu_count()): sys.exit("Cannot use more than " + str(mp.cpu_count()) + "threads")
	if (len_segment <= 0): sys.exit("Cannot use segments length that are less than 0")
	if (mapq_sr < 0): sys.exit("Cannot use MAPQ less than 0 for short reads")
	if (mapq_sr > 60): sys.exit("Cannot use MAPQ greater than 60 for short reads")
	if (mapq_lr < 0): sys.exit("Cannot use MAPQ less than 0 for long reads")
	if (mapq_lr > 60): sys.exit("Cannot use MAPQ greater than 60 for long reads")
	if (qs_sr < 0): sys.exit("Cannot use quality score less than 0 for short reads")
	if (qs_sr > 40): sys.exit("Cannot use quality score greater than 40 for short reads")
	if (qs_unmap_sr < 0): sys.exit("Cannot use quality score less than 0 for short reads")
	if (qs_unmap_sr > 40): sys.exit("Cannot use quality score greater than 40 for short reads")

	# Check that all BAM files use the same reference file
	ref_bin, ref_unbin = checkReferenceCompatibility(sr_filenames, lr_filenames, force_inter_ref)

	if (len(ref_bin.keys()) == 0): sys.exit("Input BAM files do not share any reference contigs.")

	# Open BAM files
	lr_bams = []
	sr_bams = []

	for lr_filename in lr_filenames:

		lr_bamf = pysam.AlignmentFile(lr_filename, "rb")
		lr_bams.append(lr_bamf)

	for sr_filename in sr_filenames:

		sr_bamf = pysam.AlignmentFile(sr_filename, "rb")
		sr_bams.append(sr_bamf)
	
	# Segments into different files SR and LR from same regions with (MAPQ >= mapq)
	if (nb_threads == 1):

		for chr_name, chr_len in ref_bin.items():
			
			# Segment SR and LR with MAPQ>mapq into bins
			for i in range(0, chr_len, len_segment): printBinInfo(segmentBAM(out_prefix_filename, chr_name, i, min(i + len_segment - 1, chr_len), chr_len, mapq_sr, mapq_lr, qs_sr, False))

		printBinInfo(outputUnmappedLowQualLR(out_prefix_filename, mapq_lr, ref_unbin, False)) # Output to file all LR which are either unmapped or with (MAPQ < mapq)
		printBinInfo(outputUnmappedLowQualSR(out_prefix_filename, mapq_sr, qs_unmap_sr, qs_sr, False)) # Output to file all SR which are unmapped

	else:
		
		thread_pool = mp.Pool(nb_threads) # Create pool of threads

		thread_pool.apply_async(outputUnmappedLowQualLR, args=(out_prefix_filename, mapq_lr, ref_unbin, True), callback=printBinInfo)
		thread_pool.apply_async(outputUnmappedLowQualSR, args=(out_prefix_filename, mapq_sr, qs_unmap_sr, qs_sr, True), callback=printBinInfo)

		for chr_name, chr_len in ref_bin.items():
			
			for i in range(0, chr_len, len_segment): # Segment SR and LR with MAPQ>mapq into bins

				thread_pool.apply_async(segmentBAM, args=(out_prefix_filename, chr_name, i, min(i + len_segment - 1, chr_len), chr_len, mapq_sr, mapq_lr, qs_sr, True), callback=printBinInfo)

		thread_pool.close() 
		thread_pool.join()

	# Close BAM files
	for lr_bamf in lr_bams: lr_bamf.close()
	for sr_bamf in sr_bams: sr_bamf.close()
