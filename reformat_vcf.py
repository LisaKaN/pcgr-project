import vcfpy
import argparse

def main():
    # parse input arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_vcf_file", 
        type=str, 
        help="Pass path to a vcf file"
        )
    parser.add_argument(
        "output_vcf_file", 
        type=str, 
        help="Pass path to a output file"
        )
    args = parser.parse_args()

    # read input vcf using vcfpy
    reader = vcfpy.Reader.from_path(args.input_vcf_file)

    # copy vcf header and delete header lines
    h = reader.header.copy()
    h.lines=[]

    # copy all header lines except info lines
    for l in reader.header.lines:
        if (type(l) is not vcfpy.InfoHeaderLine):
            h.add_line(l)

    # add info header lines for VAF and DP
    h.add_info_line(vcfpy.OrderedDict([('ID','TVAF'), ('Number', '.'), ('Type', 'Float'), ('Description', "Allelic fraction of alternative allele in tumor")]))
    h.add_info_line(vcfpy.OrderedDict([('ID','TDP'), ('Number', '.'), ('Type', 'Integer'), ('Description', "Read depth across variant site in tumor")]))
    
    # create new vcf file
    writer = vcfpy.Writer.from_path(args.output_vcf_file, h)

    # add all records from original vcf file to new vcf file
    for record in reader:
        dp = sum(c.data.get('DP', 0) for c in record.calls)
        vaf = list(c.data.get('AF', 0) for c in record.calls)[0][0]
        # overwrite original INFO fields and add only DP and VAF 
        record.INFO = {}
        record.INFO['TDP'] = [str(dp)]
        record.INFO['TVAF'] = [str(vaf)]
        writer.write_record(record)

if __name__ == '__main__':
    main()