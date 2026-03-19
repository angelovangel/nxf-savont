#!/usr/bin/env python3
import sys
import os
import csv

def main():
    if len(sys.argv) < 3:
        print("Usage: combine_lineage.py <output_file> <input_files...>")
        sys.exit(1)

    output_file = sys.argv[1]
    input_files = sys.argv[2:]

    tax_cols = ['tax_id', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'clade', 'superkingdom']
    combined_data = {} # (tax_tuple) -> {sample: abundance}
    samples = []

    for f in input_files:
        sample_name = os.path.basename(f).replace('_rel-abundance_lineage.tsv', '').replace('_rel-abundance.tsv', '').replace('_lineage.tsv', '')
        if sample_name not in samples:
            samples.append(sample_name)
        
        with open(f, 'r') as tsv_file:
            reader = csv.reader(tsv_file, delimiter='\t')
            for row in reader:
                if not row or row[0] == 'tax_id' or not row[0]:
                    continue
                tax_id = row[0]
                abundance = row[1] if len(row) > 1 else '0'
                if abundance == '': 
                    abundance = '0'
                
                # Default "Unknown"
                row_tax = {col: 'Unknown' for col in tax_cols}
                row_tax['tax_id'] = tax_id
                
                if len(row) >= 4:
                    names = row[2].split(';') if row[2] else []
                    ranks = row[3].split(';') if row[3] else []
                    for n, r in zip(names, ranks):
                        if not n or not r:
                            continue
                        if r == 'domain':
                            row_tax['superkingdom'] = n
                        elif r in tax_cols:
                            row_tax[r] = n
                
                tax_key = tuple(row_tax[col] for col in tax_cols)
                if tax_key not in combined_data:
                    combined_data[tax_key] = {}
                combined_data[tax_key][sample_name] = abundance

    # Sort samples for consistency
    samples.sort()

    with open(output_file, 'w', newline='') as out_file:
        writer = csv.writer(out_file, delimiter='\t')
        # Header
        writer.writerow(tax_cols + samples)
        
        # Rows
        for tax_key, sample_values in combined_data.items():
            row = list(tax_key)
            for s in samples:
                row.append(sample_values.get(s, '0'))
            writer.writerow(row)

if __name__ == "__main__":
    main()
