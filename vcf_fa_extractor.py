#!/usr/bin/env python
"""This script will extract snps from a vcf file and produce a multi-fasta
output with extracted reference and all variant sequences.
Author: Oleksandr Moskalenko <om@hpc.ufl.edu>
Version: 1.2
Date: 2014-07-25
"""
import os, sys, operator, logging

version="1.2"

def get_arguments():
    parser = argparse.ArgumentParser(usage='%(prog)s [options] -i input_file [-o output_file]', epilog="You must at least provide the input file name")
    parser.add_argument('-i', '--input', dest='infile', help="Input vcf file")
    parser.add_argument('-o', '--output', dest='outfile', help="Output fasta file")
    parser.add_argument('-t', '--table', dest='table', help="Output a tab-delimited sample data table to a file")
    parser.add_argument('-q', '--quality', dest='quality', type=int, default=-1, help="QUAL cutoff, optionsl")
    parser.add_argument('-d', '--distance', dest='distance', type=int,
            default=-1, help="Distance filter (minimal distance between bases)")
    parser.add_argument('-s', '--snpeff', dest='snpeff', help="Parse SnpEFF produced extended INFO field and output a codon alignment into the specified file.")
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
            help="Verbose output")
    args = parser.parse_args()
    if (not args.infile):
        parser.print_help()
        sys.exit(1)
    else:
        if not os.access(args.infile, os.R_OK):
            sys.exit("Cannot access input file")
    if not args.outfile:
        args.outfile = os.path.splitext(args.infile)[0] + ".fa"
        if args.verbose:
            logger.debug("It looks like you provided the input file as an argument. The output filename will be automatically generated. No quality filtering will be performed.")
            logger.debug("Input file: {0}".format(args.infile))
            logger.debug("The output file will be called %s\n" % args.outfile)
    if args.verbose:
        print "Input file: \t%s" % args.infile
        print "Output file:\t%s" % args.outfile
        if args.quality != -1:
            print "Quality cutoff: %d" % args.quality
        else:
            print "Quality cutoff is not set"
        if args.distance != -1:
            print "Distance filter: %d" % args.distance
        else:
            print "Distance filter is not set"
    return (args)

def get_header(verbose, fh):
    for line in fh:
        if line.strip().startswith('##'):
            pass
        elif line.strip().startswith('#CHROM'):
            header = line.strip()
            header_list = header.split('\t')
            header_list[0] = "CHROM"
            header_str = ", ".join(header_list)
            return header_list
        else:
            logger.error("Could not find the #CHROM header. Is this a VCF file?")
            sys.exit(1)

def parse_snpeff_info(args, data):
    # INFO=<ID=EFF,Number=.,Type=String,Description="Predicted effects for this variant.Format: 'Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_change| Amino_Acid_length | Gene_Name | Gene_BioType | Coding | Transcript | Exon  | GenotypeNum [ | ERRORS | WARNINGS ] )' ">
#    INFO AB=0;ABP=0;AC=2;AF=0.0322581;AN=62;AO=6;BasesToClosestVariant=32;CIGAR=1X;DP=10291;DPRA=1.19446;EPP=3.0103;EPPR=130.21;HWE=-0;LEN=1;MEANALT=1;MQM=32.5;MQMR=41.7045;NS=62;NUMALT=1;ODDS=5.116;PAIRED=0.333333;PAIREDR=0.769163;RO=10280;RPP=4.45795;RPPR=22031.4;RUN=1;SAP=3.0103;SRP=127.6;TYPE=snp;XAI=0;XAM=0.0188826;XAS=0.0188826;XRI=0.000124714;XRM=0.00114359;XRS=0.00101888;technology.ILLUMINA=1;EFF=SYNONYMOUS_CODING(LOW|SILENT|ggT/ggC|G91|716|Vch1786_I0103||CODING|Vch1786_I0103|1|1)
#SNPEff record:
#EFF=NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|aCc/aTc|T34I|197|Vch1786_I0077||CODING|Vch1786_I0077|1|1)
    info_list = data.split(";")
    snpeff_output = []
    table_output = []
    empty_snpeff_info_list = [''] * 13
    alt_allele = ""
    #output for the sample table file
    #Effect (Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_change| Amino_Acid_length | Gene_Name | Gene_BioType | Coding | Transcript | Exon)
    for info_item in info_list:
        if info_item.startswith("EFF"):
            if args.verbose:
                logger.debug("EFF record:")
                print info_item
            effect = info_item.split("=")[1].split("(")[0]
            effect_data_src = info_item.split("(")[1][:-1]
            effect_data_list = effect_data_src.split("|")
            allele_change_str_src = effect_data_list[2]
            [ ref_codon, alt_codon ] = allele_change_str_src.split("/")
            ref_codon.strip()
            alt_codon.strip()
            # table output - 
            # Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_change| Amino_Acid_length | Gene_Name | Gene_BioType | Coding | Transcript | Exon | SampleOne | SampleTwo
            table_output.append(effect)
            table_output.extend(effect_data_list[:10])
            if effect_data_list[10] == '1':
                table_output.extend(['yes','no'])
            elif effect_data_list[10] == '2':
                table_output.extend(['no','yes'])
            if args.verbose:
                logger.debug("Table data")
                print table_output
            break
    else:
        [ ref_codon, alt_codon ] = ['', '']
        table_output=empty_snpeff_info_list
#    print "EFF Parser output: '{0}' - '{1}' - '{2}'".format(ref_codon, alt_codon, table_output)
    return (ref_codon, alt_codon, table_output)

def parse_input(args, header, fh):
    """
    Example header:

    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
    WCH0162_SM      SJN1 1_SM       SJN15_SM        SJN12_SM        SJN01_SM
    SJN13_SM        WCH0155_SM      SJN02_SM        WCH0172_SM      SJN0 3_SM
    WCH0170_SM      SJN04_SM        WCH0149_SM      SJN05_SM        WCH0164_SM
    SJN06_SM        WCH0182_SM       SJN08_SM       SJN07_SM        WCH0151_SM
    WCH0150_SM      WCH0127_SM      SJN16_SM        SJN14_SM        SJN10_S M
    SJN09_SM        SJN17_SM

    Example INFO with the SnpEff extensions:
        INFO=<ID=EFF,Number=.,Type=String,Description="Predicted effects for this variant.Format: 'Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_change| Amino_Acid_length | Gene_Name | Gene_BioType | Coding | Transcript | Exon  | GenotypeNum [ | ERRORS | WARNINGS ] )' ">
    """
#    example = """gi|87125858|gb| CP000255.1|	2863880	.	A	G	2960.84	.	AB=0;ABP=0;AC=6;AF=0.222222;AN=27;AO=102 ;CIGAR=1X;DP=7273;DPRA=0.509969;EPP=215.9;EPPR=4.53721;HWE=-0;LEN=1;MEANALT= 1;MQM=13.8725;MQMR=54.9476;NS=27;NUMALT=1;ODDS=74.4167;PAIRED=0.990196;PAIRE DR=0.992468;RO=7169;RPP=215.9;RPPR=15070.6;RUN=1;SAP=224.5;SRP=3.11965;TYPE= snp;XAI=0.00886278;XAM=0.00929294;XAS=0.000430161;XRI=1.23666e-05;XRM=0.0008 01506;XRS=0.000789139;technology.ILLUMINA=1;BVAR	GT:GQ:DP:RO:QR:AO:QA:GL	0:5 0000:414:414:15906:0:0:0,-1431.92	0:50000:334:334:12912:0:0:0,-1162.47	1:500 00:29:0:0:29:1056:-95.4041,0	0:50000:269:269:10281:0:0:0,-925.672	1:50000:14 :0:0:14:519:-47.0807,0	0:50000:441:441:16891:0:0:0,-1520.57	0:50000:311:311: 11916:0:0:0,-1072.82	1:50000:15:1:40:14:493:-44.7221,-4	0:50000:334:333:1278 3:1:24:-2.4,-1150.85	1:50000:14:0:0:14:521:-47.2621,0	0:50000:345:345:13321: 0:0:0,-1199.28	0:50000:181:181:6920:0:0:0,-623.182	0:50000:336:336:12882:0:0 :0,-1159.76	0:50000:423:423:16229:0:0:0,-1460.99	0:50000:283:283:10902:0:0:0 ,-981.565	1:50000:16:2:65:14:529:-47.9879,-6.175	0:50000:514:514:19905:0:0:0 ,-1791.84	0:50000:377:374:14472:3:118:-11.0133,-1302.87	0:50000:368:367:1410 1:1:39:-3.9,-1269.47	0:50000:267:267:10237:0:0:0,-921.713	0:50000:226:226:87 00:0:0:0,-783.385	1:50000:10:0:0:10:362:-32.942,0	0:50000:377:376:14405:1:9: -0.9,-1296.83	0:50000:334:333:12842:1:41:-4.1,-1156.17	0:50000:262:262:10073 :0:0:0,-906.954	0:50000:495:493:19004:0:0:-3.895,-1714.43	0:50000:284:284:10 969:0:0:0,-987.596"""
    all_data = []
    all_data_table = []
    no_eff_data_sample_names = []
    output_header = list(operator.itemgetter(0,1,3,4,5)(header))
    output_header.extend(['EFF_REF', 'EFF_ALT'])
    samples = operator.itemgetter(slice(9,None))(header)
    output_header.extend(samples)
    if args.verbose:
        logger.debug("Output header")
        print ", ".join(output_header)
    if args.table:
        data_table_header = [ 'Chrom', 'Position', 'Effect', 'Effect_Impact',
            'Functional_Class', 'Codon_Change', 'Amino_Acid_change',
            'Amino_Acid_length', 'Gene_Name', 'Gene_BioType', 'Coding',
            'Transcript', 'Exon', 'SampleOne (yes/no)', 'SampleTwo (yes/no)']
        all_data_table.append(data_table_header)
        if args.verbose:
            logger.debug("Output table header")
            logger.debug(", ".join(data_table_header))
    else:
        if args.verbose:
            logger.debug("Not producing the data table")
    all_data.append(output_header)
    for line in fh:
        if not line.startswith('#'):
            try:
                all_samples_data = []
                data_table = []
                raw_record = line.strip().split('\t')
                if args.verbose:
                    logger.debug("Raw record")
                    print raw_record
                common_sample_data = list(operator.itemgetter(0,1,3,4,5)(raw_record))
                all_samples_data.extend(common_sample_data)
                data_table.extend(common_sample_data[:2])
                all_samples_data_raw = operator.itemgetter(slice(9,None))(raw_record)
                info_data_raw_str = raw_record[7]
                if args.snpeff:
                    ref_codon, alt_codon, table_output = parse_snpeff_info(args, info_data_raw_str)
                    if not ref_codon and not alt_codon:
                        no_eff_data_sample_names.append("_".join(all_samples_data[:2]))
                    data_table.extend(table_output)
                #Replace alt with SNPEFF codon in common data if asked for at index [3]
                    all_samples_data.extend([ref_codon, alt_codon])
                else:
                    all_samples_data.extend(['', ''])
                for i in all_samples_data_raw:
                    all_samples_data.append(i.strip().split(':')[0])
#                all_samples_data.append(info_data_raw_str)
                all_data.append(all_samples_data)
                all_data_table.append(data_table)
            except IOError:
                continue
            except IndexError:
                continue
    fh.close()
    if args.verbose:
        num_all_entries = len(all_data) - 1
        print "Read %d entries from the input file.\n" % num_all_entries
    return all_data, all_data_table, no_eff_data_sample_names

def write_fasta_file(verbose, filename, data):
    if verbose:
        print "Writing fasta data to {0} file".format(filename)
    try:
        fh = open(filename, 'w')
    except:
        sys.exit("Cannot open fasta alignment file for writing.")
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

def write_codon_alignment_output(verbose, filename, data):
    if verbose:
        print "Writing the fasta formatted codon alignment file to '%s' file\n" % filename
    try:
        fh = open(filename, 'w')
    except:
        sys.exit("Cannot open codon alignment file for writing.")
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

def write_data_table(args, data):
    # Chromosome | Position | Effect | Codon_Change | Amino_acid_change | Gene_Name | Sample 1 (yes/no) | Sample 2 (yes,no) 
    outfile = args.table
    try:
        fh = open(outfile, 'w')
    except:
        sys.exit("Could not open the header table file for writing.\n")
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

def convert_to_codon_alignment_fasta(verbose, specimen, filtered_data):
    samples = {}
    seq_data = {}
    for sample in specimen:
        seq_data[sample] = []
    seq_data["reference"] = []
    for entry in filtered_data:
        ref_seq = entry[5].strip()
        alt_seq = entry[6].strip()
#        print "DEBUG: EFF codons: '{0}'/'{1}'".format(ref_seq, alt_seq)
        snps = entry[7:]
        if len(specimen) != len(snps):
            print "Error: Number of specimen in the header and the sequence data for specimen {0} does not match for the following record:"
            print ", ".join(entry)
            sys.exit()
        if not ref_seq and not alt_seq:
            if verbose:
                print "Empty EFF field for:"
                print ", ".join(entry)
            continue
        seq_data["reference"].append(ref_seq)
        for sample in specimen:
            sample_id = specimen.index(sample)
            snp = snps[sample_id]
            if snp == '.':
                snp_value = ''
            else:
                try:
                    snp = int(snp)
                except:
                    sys.exit("Error: SNP call is not a '.' or an integer in the (0-n) range.")
                if snp == 0:
                    snp_value = ref_seq
                else:
                    snp_value = alt_seq
            seq_data[sample].append(snp_value)
    return seq_data

def convert_to_fasta(verbose, specimen, filtered_data):
    fasta_data = []
    seq_data = {}
    for sample in specimen:
        seq_data[sample] = []
    seq_data["reference"] = []
    for entry in filtered_data:
        ref_seq = entry[2].strip()
        alt_seq = entry[3].strip().split(',')
        snps = entry[7:]
        if len(specimen) != len(snps):
            sys.exit("Error: Number of specimen in the header and the sequence data does not match.")
        seq_data["reference"].append(ref_seq)
        for sample in specimen:
            sample_id = specimen.index(sample)
            snp = snps[sample_id]
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

def filter_data(args, source_data):
    quality = args.quality
    distance = args.distance
    reference = []
    header = source_data[0]
    original_data = source_data[1:]
    specimen = header[7:]
    sample_names = ", ".join(specimen)
    if args.verbose:
        print "Common data: %s\n" % ", ".join(header[:7])
        print "Sample names: %s\n" % sample_names
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

if __name__=='__main__':
    import sys
    logger = logging.getLogger(__name__)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s:%(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.setLevel('DEBUG')
    if sys.version_info < (2,7,0):
        print ("You need python 2.7 or later to run this script.")
        sys.exit(1)
    import argparse
    args = get_arguments()
    try:
        input_fh = open(args.infile, 'r')
    except IOError:
        print 'Cannot open input file', arg
    header = get_header(args.verbose, input_fh)
    # CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLES*
    input_fh.seek(0, 0)
    source_data, data_table, no_eff_data_samples = parse_input(args, header, input_fh)
#    if args.verbose:
#        logger.debug("All data before filtering")
#        for snp in source_data:
#            print snp
#        print "These {0} samples do not have EFF data:".format(len(no_eff_data_samples))
#        print ", ".join(no_eff_data_samples)
#        logger.debug("Data table")
#        print data_table
    if args.table:
        write_data_table(args, data_table)
    specimen, filtered_data = filter_data(args, source_data)
#    if args.verbose:
#        logger.debug('Specimen:')
#        print specimen
#        logger.debug('Filtered data:')
#        print filtered_data
#        sys.exit()
    #VERIFIED ---------------------------------------------------------
    fasta_data = convert_to_fasta(args.verbose, specimen, filtered_data)
    write_fasta_file(args.verbose, args.outfile, fasta_data)
    if args.snpeff:
        codon_alignment = convert_to_codon_alignment_fasta(args.verbose, specimen, filtered_data)
        write_codon_alignment_output(args.verbose, args.snpeff, codon_alignment)
    if args.verbose:
        print "Done processing the vcf file {0}. Good bye!\n".format(args.infile)
    sys.exit(0)
