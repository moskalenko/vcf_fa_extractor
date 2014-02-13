#!/usr/bin/env python
"""This script will extract snps from a vcf file and produce a multi-fasta
output with extracted reference and all variant sequences.
Author: Oleksandr Moskalenko <om@hpc.ufl.edu>
Version: 1.1
Date: 2014-02-13
"""

import os, sys, operator

def get_arguments():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('-i', '--input', dest='infile', help="Input vcf file", metavar="INFILE")
    parser.add_option('-o', '--output', dest='outfile', help="Output fasta file", metavar="OUTFILE")
    parser.add_option('-q', '--quality', dest='quality', type="int", help="QUAL cutoff, optionsl", metavar="QUALITY")
    parser.add_option('-d', '--distance', dest='distance', type="int", help="Distance filter (minimal distance between bases)", metavar="DISTANCE")
    parser.add_option('-v', '--verbose', action='store_true', default=False, help="Verbose output", metavar="VERBOSE")
    (opts, args) = parser.parse_args()
    infile_arg = ''
    outfile_arg = ''
    verbose_arg = False
    if opts.verbose:
        verbose_arg = opts.verbose
    if (not opts.infile) and (len(args) == 0):
        print "You must provide at least the input file name.\n"
        usage()
        sys.exit("Wrong input(s) detected\n")
    if len(args) == 1:
        if verbose_arg:
            print "It looks like you provided the input file as an argument. The output filename will be automatically generated. No quality filtering will be performed.\n"
            print "Input file: %s.\n" % args[0]
        infile_arg = args[0]
        outfile_arg = os.path.splitext(infile_arg)[0] + ".fa"
        if verbose_arg:
            print "The output file will be called %s\n" % outfile_arg
    if not infile_arg:
        if opts.infile:
            infile_arg = opts.infile
        else:
            usage()
            sys.exit("Please provide an input file name\n")
    if opts.outfile:
        outfile_arg = opts.outfile
    quality_arg = -1
    distance_arg = -1
    if not outfile_arg:
        outfile_arg = os.path.splitext(infile_arg)[0] + ".fa"
    if opts.quality:
        quality_arg = opts.quality
    if opts.distance:
        distance_arg = opts.distance
    if verbose_arg:
        print "Input file: %s" % infile_arg
        print "Output file: %s" % outfile_arg
        if quality_arg != -1:
            print "Quality cutoff: %d" % quality_arg
        else:
            print "Quality cutoff is not set"
        if distance_arg != -1:
            print "Distance filter: %d" % distance_arg
        else:
            print "Distance filter is not set"
    return (infile_arg, outfile_arg, quality_arg, distance_arg, verbose_arg)

def usage():
    print """Usage: ./vcf_extract_taj.py [options]

                Options: 
                    -h, --help            show this help message and exit
                    -i INFILE, --input=INFILE
                                        Input vcf file
                    -o OUTFILE, --output=OUTFILE
                                        Output fasta file

                Alternatively provide input and output file names as arguments:
                     ./vcf_extract_taj.py input.vcf output.fa
          """

def get_header(verbose, infile):
    fh = open(infile, 'r')
    for line in fh:
        if line.strip().startswith('##'):
            pass
        elif line.strip().startswith('#CHROM'):
            header = line.strip()
            return (header, fh)
        else:
            return ('', fh)

def parse_input(verbose, header, fh):
    """
    Example header:

    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
    WCH0162_SM      SJN1 1_SM       SJN15_SM        SJN12_SM        SJN01_SM
    SJN13_SM        WCH0155_SM      SJN02_SM        WCH0172_SM      SJN0 3_SM
    WCH0170_SM      SJN04_SM        WCH0149_SM      SJN05_SM        WCH0164_SM
    SJN06_SM        WCH0182_SM       SJN08_SM       SJN07_SM        WCH0151_SM
    WCH0150_SM      WCH0127_SM      SJN16_SM        SJN14_SM        SJN10_S M
    SJN09_SM        SJN17_SM

    Example data:
    """
#    example = """gi|87125858|gb| CP000255.1|	2863880	.	A	G	2960.84	.	AB=0;ABP=0;AC=6;AF=0.222222;AN=27;AO=102 ;CIGAR=1X;DP=7273;DPRA=0.509969;EPP=215.9;EPPR=4.53721;HWE=-0;LEN=1;MEANALT= 1;MQM=13.8725;MQMR=54.9476;NS=27;NUMALT=1;ODDS=74.4167;PAIRED=0.990196;PAIRE DR=0.992468;RO=7169;RPP=215.9;RPPR=15070.6;RUN=1;SAP=224.5;SRP=3.11965;TYPE= snp;XAI=0.00886278;XAM=0.00929294;XAS=0.000430161;XRI=1.23666e-05;XRM=0.0008 01506;XRS=0.000789139;technology.ILLUMINA=1;BVAR	GT:GQ:DP:RO:QR:AO:QA:GL	0:5 0000:414:414:15906:0:0:0,-1431.92	0:50000:334:334:12912:0:0:0,-1162.47	1:500 00:29:0:0:29:1056:-95.4041,0	0:50000:269:269:10281:0:0:0,-925.672	1:50000:14 :0:0:14:519:-47.0807,0	0:50000:441:441:16891:0:0:0,-1520.57	0:50000:311:311: 11916:0:0:0,-1072.82	1:50000:15:1:40:14:493:-44.7221,-4	0:50000:334:333:1278 3:1:24:-2.4,-1150.85	1:50000:14:0:0:14:521:-47.2621,0	0:50000:345:345:13321: 0:0:0,-1199.28	0:50000:181:181:6920:0:0:0,-623.182	0:50000:336:336:12882:0:0 :0,-1159.76	0:50000:423:423:16229:0:0:0,-1460.99	0:50000:283:283:10902:0:0:0 ,-981.565	1:50000:16:2:65:14:529:-47.9879,-6.175	0:50000:514:514:19905:0:0:0 ,-1791.84	0:50000:377:374:14472:3:118:-11.0133,-1302.87	0:50000:368:367:1410 1:1:39:-3.9,-1269.47	0:50000:267:267:10237:0:0:0,-921.713	0:50000:226:226:87 00:0:0:0,-783.385	1:50000:10:0:0:10:362:-32.942,0	0:50000:377:376:14405:1:9: -0.9,-1296.83	0:50000:334:333:12842:1:41:-4.1,-1156.17	0:50000:262:262:10073 :0:0:0,-906.954	0:50000:495:493:19004:0:0:-3.895,-1714.43	0:50000:284:284:10 969:0:0:0,-987.596"""
    all_data = []
    header_list = header.split('\t')
    header_list[0] = "CHROM"
    short_header = list(operator.itemgetter(0,1,3,4,5)(header_list))
    samples = operator.itemgetter(slice(9,None))(header_list)
    short_header.extend(samples)
    all_data.append(short_header)
    for line in fh:
        try:
            raw_input_data = line.strip().split('\t')
            header_data = list(operator.itemgetter(0,1,3,4,5)(raw_input_data))
            sample_data = []
            raw_sample_data = operator.itemgetter(slice(9,None))(raw_input_data)
            for i in raw_sample_data:
                sample_data.append(i.strip().split(':')[0])
            header_data.extend(sample_data)
            all_data.append(header_data)
        except IOError:
            continue
        except IndexError:
            continue
    fh.close()
    if verbose:
        num_all_entries = len(all_data) - 1
        print "Read %d entries from the input file.\n" % num_all_entries
    return all_data

def write_output(verbose, outfile, data):
    fh = open(outfile, 'w')
    for line in data:
        fh.write(",".join(line) + os.linesep)
    fh.close()

def quality_filter(verbose, quality, original_data):
    qual = float(quality)
    num_original_entries = len(original_data)
    if verbose:
        print "Received %d entries for quality filtering with a %.1f cutoff.\n" % (num_original_entries, qual)
    filtered_data = []
    if verbose:
        print "Quality cutoff: %.1f\n" % qual
    for entry in original_data:
        entry_quality = str(entry[4])
        entry_quality = entry_quality.replace(' ','')
#        print "Quality: '%s'\n" % entry_quality
        if float(entry_quality) > qual:
            filtered_data.append(entry)
#    print "\nQuality filtered data:\n"
#    print quality_filtered_data
    num_filtered_entries = len(filtered_data)
    if verbose:
        print "After the quality filtering %d entries remain.\n" % (num_filtered_entries)
    return filtered_data

def distance_filter(verbose, distance, original_data):
    num_original_entries = len(original_data)
    if verbose:
        print "Received %d entries for distance filtering with a %d interval.\n" % (num_original_entries, distance)
    filtered_data = []
    discarded_data = []
    distances = []
    discarded_distances = []
    first = original_data[0]
    second = original_data[1]
    last = original_data[-1]
    before_last = original_data[-2]
    before_before_last = original_data[-3]
    if int(second[1]) > (int(first[1]) + distance):
        filtered_data.append(first)
        distances.append(first[1])
    else:
        discarded_data.append(first)
        discarded_distances.append(first[1])
    tail = original_data[0]
    mid = original_data[1]
    head = original_data[2]
    for entry in original_data[3:]:
        if (int(tail[1]) + distance) <  int(mid[1]) < (int(head[1]) - distance):
            filtered_data.append(mid)
            distances.append(mid[1])
        else:
            discarded_data.append(mid)
            discarded_distances.append(mid[1])
#            if verbose:
#                print "Distance filter triggered on:"
#            print "Previous: ", tail
#                print "Discarded: ", mid
#            print "Next: ", head
        tail = mid
        mid = head
        head = entry
    if int(before_last[1]) > (int(before_before_last[1]) + distance):
        filtered_data.append(before_last)
        distances.append(before_last[1])
    else:
        discarded_data.append(before_last)
        discarded_distances.append(before_last[1])
    if int(last[1]) > (int(before_last[1]) + distance):
        filtered_data.append(last)
        distances.append(last[1])
    else:
        discarded_data.append(last)
        discarded_distances.append(last[1])
    num_filtered_entries = len(filtered_data)
    num_discarded_entries = len(discarded_data)
    if verbose:
        print "After the distance filtering %d entries remain. %d entries were discarded\n" % (num_filtered_entries, num_discarded_entries)
    return filtered_data

def convert_to_fasta(verbose, header, filtered_data):
    samples = {}
    seq_data = {}
    for sample in specimen:
        seq_data[sample] = []
    seq_data["reference"] = []
    for entry in filtered_data:
        ref_seq = entry[2].strip()
        alt_seq = entry[3].strip().split(',')
        snps = entry[5:]
        if len(header) != len(snps):
            sys.exit("Error: Number of specimen in the header and the sequence data does not match.")
        seq_data["reference"].append(ref_seq)
        for sample in header:
            sample_id = header.index(sample)
            snp = snps[sample_id]
#            if verbose:
#                print "SNP: %s" % snp
            if snp == '.':
                snp_value = '?'
            else:
                try:
                    snp = int(snp)
                except:
                    sys.exit("Error: SNP call is not a '.' or an integer in the (0-n) range.")
                if snp == 0:
                    snp_value = ref_seq
                else:
                    alt_snp = snp - 1
                    snp_value = alt_seq[alt_snp]
            seq_data[sample].append(snp_value)
    return seq_data

def filter_data(verbose, source_data, quality, distance):
    reference = []
    header = source_data[0]
    original_data = source_data[1:]
    specimen = header[5:]
    sample_names = ", ".join(specimen)
    if verbose:
        print "Sample names: %s\n" % sample_names
#        echo = 0
#        print "Samples (1-5):"
#        while True:
#            if echo < 5:
#                print original_data[echo]
#                echo += 1
#            else:
#                break
    #Quality filter
    if quality != -1:
        quality_filtered_data = quality_filter(verbose, quality, original_data)
    else:
        quality_filtered_data = original_data
    #Distance filter
    if distance != -1:
        distance_filtered_data = distance_filter(verbose, distance, quality_filtered_data)
    else:
        distance_filtered_data = quality_filtered_data
    return (specimen, distance_filtered_data)

def write_output(verbose, filename, data):
    if verbose:
        print "Writing the fasta formatted data to the '%s' file\n" % filename
    try:
        fh = open(filename, 'w')
    except:
        sys.exit("Cannot open output file for writing.")
    reference_list = data["reference"]
    reference_fa = "".join(reference_list)
    fh.write('>reference\n')
    fh.write(reference_fa)
    fh.write("\n")
    for specimen in data:
        if specimen != 'reference':
            specimen_name = ">" + specimen + "\n"
            specimen_fa = "".join(data[specimen]) + "\n"
            fh.write(specimen_name)
            fh.write(specimen_fa)
    fh.close()

if __name__=='__main__':
    infile, outfile, quality, distance, verbose = get_arguments()
    if verbose:
        print "Verbose output has been requested.\n"
#    print "Input: '%s'\nOutput: '%s'\n" % (infile_arg, outfile_arg)
    (header, data_handle) = get_header(verbose, infile)
    source_data = parse_input(verbose, header, data_handle)
    specimen, filtered_data = filter_data(verbose, source_data, quality, distance)
    fasta_data = convert_to_fasta(verbose, specimen, filtered_data)
    write_output(verbose, outfile, fasta_data)
    if verbose:
        print "Done processing the vcf file. Good bye!\n"
