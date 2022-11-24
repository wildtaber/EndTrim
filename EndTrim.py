#!/usr/bin/python3.8
# _*_ coding:utf-8 _*_
"""
Update log for v2:
1. Trim 3' end umi adaptor with 12bp sequence in another read's 5' end.
2. Use regex to find mismatched sequence.

"""
import datetime
import gzip
import sys
from argparse import ArgumentParser
import re
import Levenshtein

umi_l = ['GGTTACT', 'CATAGGA', 'GAATCTG', 'AATCCGA', 'GTAGTAT', 'CAGTAGT', 'TAACCGT', 'CGTAATC', 'TATAGCG', 'AGCTAAC',
        'GACGTAA', 'AAGAGTC', 'TGTAAGA', 'GTTACAA', 'GAGTACA', 'CGGTTAT', 'CTTCAAC', 'CCATTAG', 'AAGGTTG', 'CGATATA',
        'ACATAGA', 'GAGTTAG', 'AACGCAT', 'ACAATGT', 'GGATAAG', 'AGATTGC', 'TGAACTT', 'ACTTATC', 'TAGATGC', 'TCTGCTT',
        'ATCGTCT', 'GTTAACG', 'TCGTATA', 'TGATAGT', 'TTCTTCC', 'ACAACAC', 'CATGTAG', 'TGCATGT', 'GAATGGA', 'TGACCAA',
        'CCTAGTT', 'TGTTCAG', 'AGATCCA', 'TATGGAA', 'CTTGGAT', 'AAGACAG', 'ATCAAGC', 'GAGACTT', 'TCTCGTA', 'TCATGTG',
        'ATCTCAG', 'CTGTATG', 'ACTCCAT', 'CTAAGTC', 'AACAGGT', 'CTTCTTG', 'CTATAGC', 'CTACCTA', 'GTGAGTA', 'TCGTTAC',
        'GTCATGA', 'CGAATAA', 'GTGTAAC', 'TAAGGTC', 'AGACACT', 'AACCGAA', 'AGTTAGG', 'CACTATC', 'CTGATAC', 'TGCCATA',
        'AACACCA', 'ACATGCT', 'CTCCTAA', 'TGAATCG', 'TCGAAGT', 'CATTACG', 'ACCGATA', 'GATGCTA', 'GTTGAGA', 'AGTAGCA',
        'TTGGACA', 'GTAAGCT', 'TATCAGC', 'TAACGCA', 'TCTCTAG', 'AATCGCT', 'GTACAGT', 'CTTATCA', 'TCTATCT', 'GTTCCTT',
        'CTATTCT', 'TTCGCTA', 'TCTGTGA', 'TAAGTCT', 'TTCAACT', 'AAGTCTA', 'GTTATTC', 'GAAGATT', 'GTCTATT', 'GATAGAC']


def gzfq_iterator(read1_fastq, read2_fastq):
    read1_readline = read1_fastq.readline
    read2_readline = read2_fastq.readline

    while True:
        read1_line = bytes.decode(read1_readline())
        read2_line = bytes.decode(read2_readline())

        if not read1_line and read2_line:
            return
        if read1_line[0] == '@' and read2_line[0] == '@':
            break
        if isinstance(read1_line[0], int) or isinstance(read2_line[0], int):
            raise ValueError("FASTQ files may contain binary information or are compressed")

    while read1_line and read2_line:

        if read1_line[0] != '@' or read2_line[0] != '@':
            print(read1_line, read2_line)
            raise ValueError("Records in FASTQ files should start with a '@' character. \
                                Files may be malformed or out of synch.")

        title_read1_line = read1_line[1:].rstrip()
        title_read2_line = read2_line[1:].rstrip()

        read1_seq_string = bytes.decode(read1_readline()).rstrip()
        read2_seq_string = bytes.decode(read2_readline()).rstrip()

        while True:
            read1_line = bytes.decode(read1_readline())
            read2_line = bytes.decode(read2_readline())

            if not read1_line and read2_line:
                raise ValueError("End of file without quality information. Files may be malformed or out of synch")
            if read1_line[0] == '+' and read2_line[0] == '+':
                break

            read1_seq_string += read1_line.rstrip()
            read2_seq_string += read2_line.rstrip()

        read1_quality_string = bytes.decode(read1_readline()).rstrip()
        read2_quality_string = bytes.decode(read2_readline()).rstrip()

        while True:
            read1_line = bytes.decode(read1_readline())
            read2_line = bytes.decode(read2_readline())

            if not read1_line or read2_line:
                break  # end of file
            if (read1_line[0] == '@' and read2_line[0] == '@' and read1_line.isalpha() is not True
                    and read2_line.isalpha() is not True):
                break

            read1_quality_string += read1_line.rstrip()
            read2_quality_string += read2_line.rstrip()

        yield (title_read1_line, title_read2_line,
               read1_seq_string, read2_seq_string,
               read1_quality_string, read2_quality_string)

        # raise StopIteration


def rename_title(read_title, read_tag):
    # This function renames the header with the formatting of
    # *header coordinates,etc*, *tag from read1*, *tag from read2*,
    # *read designation from original header (for paired reads)*
    title_l = read_title.split(" ")
    illumina = title_l[0].split(':')
    if len(illumina) == 7:
        # Illumina CASAVA >=1.8
        # e.g. @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACGCC
        return "%s_%s %s" % (title_l[0], read_tag, title_l[1])
    elif len(illumina) == 5:
        # Illumina CASAVA >=1.4?
        # e.g. @HWUSI-EAS100R:6:73:941:1973#ATCGAT/1
        title = read_title.split('#')
        return "%s_%s#%s" % (title[0], read_tag, title[1])
    else:
        raise ValueError("Unknown read name format: %s" % read_title)


def completed(sequence):
    base_dic = {'A': "T", "T": "A", 'C': "G", "G": "C", "N": "N",
                'a': "t", "t": "a", 'c': "g", "g": "c", "n": "n"}
    rc_seq = ''
    for i in sequence.strip():
        rc_seq += base_dic[i]
    return rc_seq


def get_regex(s):
    regex_s = ''
    for i in range(len(s)):
        str_list = list(s)
        str_list[i] = "."
        if regex_s == "":
            regex_s = "".join(str_list)
        else:
            regex_s += "|"
            regex_s += "".join(str_list)
    return regex_s


def umi_trim(read_seq, etag):
    read_seq_r = read_seq[::-1]  # reverse seq
    etag_c = completed(etag)
    read_trim = read_seq

    etag_i = read_seq_r.find(etag_c)
    if etag_i != -1:
        read_trim = read_seq_r[etag_i:][::-1]
    else:
        for i in range(len(read_seq_r)-11):
            if Levenshtein.distance(read_seq_r[i:i+12], etag_c) <= 1:
                read_trim = read_seq_r[i:][::-1]
                break

    return read_trim


def open_fastq(infile, outfile1, outfile2):
    if infile.endswith(".gz"):
        in_fh = gzip.open(infile, 'rb')
    else:
        in_fh = open(infile, 'r')
    out_fh1 = open(outfile1, 'w')
    out_fh2 = open(outfile2, 'w')
    return in_fh, out_fh1, out_fh2


def umi_compare(umi1, umi2):
    mismatch = 0
    umi_len = min(len(umi1), len(umi2))
    for i in range(umi_len):
        if umi1[i] != umi2[i]:
            mismatch += 1
    return mismatch


def get_umi_name(umi):
    umi_name = -1
    if umi in umi_l:
        umi_name = umi
    else:
        for i in umi_l:
            if Levenshtein.distance(umi, i) <= 1:
                umi_name = i
                break
    return umi_name


def main():
    parser = ArgumentParser()
    parser.add_argument('--infile1', dest='infile1', help='Path to FASTQ file for Read 1.', required=True)
    parser.add_argument('--infile2', dest='infile2', help='Path to FASTQ file for Read 2.', required=True)
    parser.add_argument('--outprefix', dest='outfile', help='Prefix for output files. \
     Will prepend onto file name of ".fq.smi"', required=True)
    parser.add_argument('--taglen', dest='taglen', type=int, default=7,
                        help='Length in bases of the duplex tag sequence.[7]')
    parser.add_argument('--etag_len', dest='etaglen', type=int, default=12,
                         help='Length in bases of the duplex tag sequence.[12]')
    parser.add_argument('--spacerlen', dest='spclen', type=int, default=1,
                        help='Length in bases of the spacer sequence between duplex tag \
                         and the start of target DNA. [1]')
    parser.add_argument('--readout', dest='readout', type=int, default=1000000,
                        help='How many reads are processed before progress is reported. [1000000]')

    ags = parser.parse_args()
    now = datetime.datetime.now()

    sys.stderr.write("*" * 100)
    sys.stderr.write("\n%s\n3' UMI trimming start \n" % now)
    sys.stderr.write("Read1: %s\n" % ags.infile1)
    sys.stderr.write("Read2: %s\n\n" % ags.infile2)

    (read1_fastq, read1_output, badU_out1) = open_fastq(ags.infile1, ags.outfile + '_R1_t.fq', ags.outfile + '_R1_bu.fq')
    (read2_fastq, read2_output, badU_out2) = open_fastq(ags.infile2, ags.outfile + '_R2_t.fq', ags.outfile + '_R2_bu.fq')
    readctr = 0
    badumi = 0

    for read1_title, read2_title, read1_seq, read2_seq, read1_qual, read2_qual \
            in gzfq_iterator(read1_fastq, read2_fastq):
        readctr += 1
        umi_regex = '(^.{%d})(.{1})(.{%d})' % (ags.taglen, ags.etaglen)  # 7umi,N,12etag
        rule = re.compile(umi_regex)
        umi_tag1 = re.findall(rule, read1_seq)
        umi_tag2 = re.findall(rule, read2_seq)

        if umi_tag1 and umi_tag2:
            #time = datetime.datetime.now()
            umi_tag = str(get_umi_name(umi_tag1[0][0])) + str(get_umi_name(umi_tag2[0][0]))
            #print('umi_time', datetime.datetime.now() - time)
            if "-1" not in umi_tag:
                fam_tag = umi_tag
                read1_title = rename_title(read1_title, fam_tag)
                read2_title = rename_title(read2_title, fam_tag)
                #time2 = datetime.datetime.now()
                #print(read1_seq[(ags.taglen + 1):], umi_tag2[0][2])
                #print( read2_seq[(ags.taglen + 1):], umi_tag1[0][2] )
                read1_trim = umi_trim(read1_seq[(ags.taglen + 1):], umi_tag2[0][2])
                read2_trim = umi_trim(read2_seq[(ags.taglen + 1):], umi_tag1[0][2])
                #print('trim_time', datetime.datetime.now() - time2)

                read1_output.write(
                    '@%s\n%s\n+\n%s\n' % (read1_title, read1_trim,
                                          read1_qual[ags.taglen + 1:ags.taglen + 1 + len(read1_trim)]))
                read2_output.write(
                    '@%s\n%s\n+\n%s\n' % (read2_title, read2_trim,
                                          read2_qual[ags.taglen + 1:ags.taglen + 1 + len(read2_trim)]))
            else:
                badU_out1.write(
                    '@%s\n%s\n+\n%s\n' % (read1_title, read1_seq, read1_qual))
                badU_out2.write(
                    '@%s\n%s\n+\n%s\n' % (read2_title, read2_seq, read2_qual))
                badumi += 1
        else:
            badU_out1.write(
                '@%s\n%s\n+\n%s\n' % (read1_title, read1_seq, read1_qual) )
            badU_out2.write(
                '@%s\n%s\n+\n%s\n' % (read2_title, read2_seq, read2_qual) )
            badumi += 1

        if readctr % ags.readout == 0:
            sys.stderr.write("%s\n" % datetime.datetime.now())
            sys.stderr.write("Total reads processed: %d\n" % readctr)
            sys.stderr.write("Bad UMI: %d\n\n" % badumi)

    read1_fastq.close()
    read2_fastq.close()
    read1_output.close()
    read2_output.close()
    badU_out1.close()
    badU_out2.close()

    sys.stderr.write("%s\n" % datetime.datetime.now())
    sys.stderr.write("Total reads processed: %s\n" % readctr)
    sys.stderr.write("Bad UMI: %d\n" % badumi)
    sys.stderr.write("3' UMI trimming completed: %s\n" % (datetime.datetime.now() - now))
    sys.stderr.write("*" * 100 + '\n\n')


if __name__ == "__main__":
    main()
