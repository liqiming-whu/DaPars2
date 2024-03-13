import argparse


def extract_total_reads(input_flagstat_file):
    num_line = 0
    total_reads = None
    #print input_flagstat_file
    for line in open(input_flagstat_file,'r'):
        num_line += 1
        if num_line == 7:
            assert "mapped" in line
            total_reads = line.strip().split(' ')[0]
            break
    return total_reads

def extract_read_depth(input_wig_files, input_flagstat_files, output_file):
    fho = open(output_file,'w')
    for wig_file, flagstat_file in zip(input_wig_files,input_flagstat_files):
        read_depth = extract_total_reads(flagstat_file)
        fho.write(f"{wig_file}\t{read_depth}\n")
    fho.close()    


def extract_read_depth_main():
    parser = argparse.ArgumentParser(description='Extract read depth of each sample',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--wigfiles',dest="wigfiles", nargs="+", required=True, help="wig files")
    parser.add_argument('--flagstat-files',dest="flagstat", nargs="+", required=True, help="flagstat files")
    parser.add_argument('--output',dest="output", required=True, help="the final output file with read depth of each sample")
    args = parser.parse_args()

    samples = args.wigfiles
    flagstat = args.flagstat
    assert len(samples) == len(flagstat)
    
    extract_read_depth(samples, flagstat, args.output)
    

if __name__ == '__main__':
    extract_read_depth_main()
