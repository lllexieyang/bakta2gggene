# Usage:
# Extract "cds" annotation of contigs containing specific gene from bakta results (.tsv)

# EXAMPLES:
# 1. Basic usage
#   Extracting [gene_name] from all .tsv files in ./
#      python3 bakta2gggene.py -g [gene_name]
#   Extracting [gene_name] from input.tsv in ./
#       python3 bakta2gggene.py -i input.tsv -g [gene_name]
#   Extracting [gene_name] from input.tsv in ./ and save results in ./result/
#       python3 bakta2gggene.py -i input.tsv -o result -g [gene_name]
#   Extracting [gene_name] from input1.tsv and input2.tsv in ./
#       python3 bakta2gggene.py -i input1.tsv input2.tsv -g [gene_name]
#   Extracting [gene_name] from all .tsv files in [directory of bakta tsv results]
#       python3 bakta2gggene.py -d [directory of bakta tsv results] -g [gene_name]
#   Extracting [gene_name] from input.tsv files in [directory of bakta tsv results]
#       python3 bakta2gggene.py -d [directory of bakta tsv results] -i input.tsv -g [gene_name]
# 2. Merge results into one and generate gggenes input file
#   Note: if -d, merge; if -i, no merge; if -m/--merge, merge.
#       python3 bakta2gggene.py -i input1.tsv input2.tsv -g [gene_name] -m
# 3. Maximum distance from the target gene
#   Note: default 5kb
#      python3 bakta2gggene.py -i input1.tsv input2.tsv -g [gene_name] --merge --max 10000
# 4. Exact match of gene name
#    Extracting genes start with "mcr-1.1" using "-g mcr-1.1" (mcr-1.1, mcr-1.10, et al.)
#    or with exactly the name provided using "-g mcr-1.1 -e" (mcr-1.1 only)

# Copy right reserved : Lu Yang (yanglu2016@cau.edu.cn)
# Last change: Sep 11 2023
# Version 1.0

import os
import re
import csv
import sys
import argparse

def get_arguments():
    parser = argparse.ArgumentParser(description="Bakta Parser")

    parser.add_argument("-i", "--input", action="append", dest="input", help="bakta tsv file name(s)", default=[], nargs='+')
    parser.add_argument("-d", "--directory", action="store", dest="direct", help="directory of bakta tsv files", default='./')
    parser.add_argument("-g", "--gene", action="store", dest="gene", help="gene to extract", default="")
    parser.add_argument("-o", "--output", action="store", dest="output", help="output path", default="./")
    parser.add_argument("-m", "--merge", action="store_true", dest="merge", help="merge results into one", default="")
    parser.add_argument("--max", action="store", dest="max",
                        help="Maximum distance from the target gene, default is 5kb", default="5000")
    parser.add_argument("-e", "--exact", action="store_true", dest="exact",
                        help="gene name provided exactly", default=False)

    return parser.parse_args()


def find_gene_in_bakta(file, gene, out, max_dis, exact):
    output_name = os.path.join(out, f"{gene}_from_{os.path.basename(file)}")
    with open(file, 'r') as input_file, open(output_name, 'w', newline='') as output_file:
        # Skip first two lines
        next(input_file)
        next(input_file)
        reader = csv.DictReader(input_file, delimiter='\t')
        writer = csv.DictWriter(output_file, fieldnames=reader.fieldnames, delimiter='\t')
        writer.writeheader()
        temp = {}
        sequence_list = []
        for row in reader:
            sequence_id = row['#Sequence Id']
            start = int(row['Start'])
            stop = int(row['Stop'])
            max_dis_start = start - int(max_dis)
            max_dis_stop = stop + int(max_dis)

            temp.setdefault(sequence_id, []).append({'row': row, 'start': start, 'stop': stop})

            # Find sequence_id and its start/stop position of target gene
            gene_found = 0
            if exact:
                if row['Gene'] == gene:
                    gene_found = 1
                    print("exact:", exact)
            elif row['Gene'].startswith(gene):
                gene_found = 1
                print(row['Gene'])

            if gene_found:
                if any(entry['sequence_id'] == sequence_id for entry in sequence_list):
                    for entry in sequence_list:
                        if entry['sequence_id'] == sequence_id:
                            entry['max_dis_start'] = min(entry['max_dis_start'], max_dis_start)
                            entry['max_dis_stop'] = max(entry['max_dis_stop'], max_dis_stop)
                        #print(sequence_id, "acceptable range changed: from", entry['max_dis_start'], "to", entry['max_dis_stop'])
                else:
                    sequence_list.append(
                        {'sequence_id': sequence_id, 'max_dis_start': max_dis_start, 'max_dis_stop': max_dis_stop})
                    #print(sequence_id, "acceptable range: from", max_dis_start, "to", max_dis_stop)

        # Check distance of other genes in the same conitg with target gene
        contigs = [entry['sequence_id'] for entry in sequence_list]
        if sequence_list != []:
            for item in sequence_list:
                for entry in temp.get(item['sequence_id'], []):
                    if item['max_dis_start'] <= entry['start'] <= item['max_dis_stop'] and \
                            item['max_dis_start'] <= entry['stop'] <= item['max_dis_stop']:
                        writer.writerow(entry['row'])
            print(contigs, "were extracted from", os.path.basename(file), "as", os.path.basename(output_name))
        else:
            print("No contigs containing", gene, "were found in", file)
            output_file.close()
            os.remove(output_name)
            return None
    return os.path.basename(file), contigs, output_name


def merge_results(result_list, out):
    result_names = [info['result_name'] for info in result_list.values() if 'result_name' in info]
    with open(out, 'w', newline='') as output_file:
        # Extract the column names from the first result file as fieldnames
        with open(result_names[0], 'r') as first_input_file:
            first_reader = csv.DictReader(first_input_file, delimiter='\t')
            fieldnames = ['file_name'] + first_reader.fieldnames
            writer = csv.DictWriter(output_file, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()

            for result_name in result_names:
                with open(result_name, 'r') as input_file:
                    reader = csv.DictReader(input_file, delimiter='\t')
                    for row in reader:
                        # Extract the origin file name from the result_name
                        file_name = re.search(r'from_(.*?)\.tsv', result_name).group(1)
                        row['file_name'] = file_name
                        writer.writerow(row)
    return 0


def generete_gggene_input(result, out):
    with open(result, 'r') as input_file, open(out, 'w', newline='') as output_file:
        reader = csv.DictReader(input_file, delimiter='\t')
        fieldnames = ['molecule', 'start', 'end', 'orientation', 'gene', 'product']
        writer = csv.DictWriter(output_file, fieldnames=fieldnames, delimiter=',')
        writer.writeheader()

        for row in reader:
            file_name = row['file_name']
            molecule = row['file_name'].replace('.tsv', '') + '_' + row['#Sequence Id']
            start = int(row['Start'])
            end = int(row['Stop'])
            gene = row['Gene']
            orientation = row['Strand']
            product = row['Product']
            gene_type = row['Type']
            if gene_type == 'cds':
                if orientation == '+':
                    writer.writerow({'file_name': file_name, 'molecule': molecule, 'start': start, 'end': end, 'gene': gene, 'orientation': orientation, 'product': product})
                else:
                    writer.writerow({'file_name': file_name, 'molecule': molecule, 'start': end, 'end': start, 'gene': gene, 'orientation': orientation, 'product': product})


def main():
    options = get_arguments()

    # Check if target gene was provided
    if not options.gene:
        print("Error: Please specify the gene using the -g option.")
        sys.exit(1)

    # Get a list of files in the input directory
    input_dir = options.direct
    if options.input:
        files = [os.path.join(input_dir, f) for f in options.input[0]]
        merge = False
    else:
        # All .tsv files except for output files of this script
        files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith('.tsv')
                 and '_from' not in f and '_record' not in f]
        merge = True
    print("Extracting annotation results of contigs containing\033[1;31m", options.gene, "\033[0mfrom\033[1;31m", len(files), "\033[0mfiles.")

    # Check if the output directory exists, and create it if it doesn't.
    if options.output:
        if os.path.exists(options.output):
            print("Results will be saved in\033[1;31m", options.output, "\033[0m")
        else:
            os.makedirs(options.output)
            print("Create a new output directory:\033[1;31m", options.output, "\033[0m")

    # Extract records from bakta results
    records = {}
    for file_name in files:
        if os.path.exists(file_name):
            record = find_gene_in_bakta(file_name, options.gene, options.output, options.max, options.exact)
            if record:
                records[len(records) + 1] = {
                    'file_name': record[0],
                    'sequence_id': record[1],
                    'result_name': record[2]}

    # Save info of records extracted
    print("Contigs containing\033[1;31m", options.gene, "\033[0mwere found in\033[1;31m", len(records),
          "\033[0mfiles.")
    if records:
        records_name = os.path.join(options.output, f"{options.gene}_records.tsv")
        with open(records_name, 'w', newline='') as extracted:
            fieldnames = ['file_name', 'sequence_id', 'result_name']
            writer = csv.DictWriter(extracted, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            for record in list(records.values()):
                writer.writerow(record)
            print("Extracting records saved in\033[1;31m", records_name, "\033[0m")
    else:
        merge = 0
        print("Merging operation will not be performed.")

    # Merge records into one file
    if options.merge:
        merge = options.merge
    if merge:
        merged_name = os.path.join(f"{options.output}/{options.gene}_from_bakta.tsv")
        merge_results(records, merged_name)
        print("Results merged ->\033[1;31m", merged_name, "\033[0m")

        # Generate gggene input file
        gggene_out = os.path.join(f"{options.output}/{options.gene}_to_gggene.csv")
        generete_gggene_input(merged_name, gggene_out)
        print("gggene input file generated ->\033[1;31m", gggene_out, "\033[0m")


    return 0


if __name__ == '__main__':
    main()
