#!/usr/bin/env python3
import sys
import os
import json
import csv
import re

def natural_sort_key(s):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split('([0-9]+)', str(s))]

def parse_rel_abundance(file_path):
    data = []
    try:
        with open(file_path, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                if row.get('tax_id') in ['unmapped', 'mapped_filtered', 'mapped_unclassified']:
                    continue
                abundance = float(row.get('abundance', 0))
                if abundance > 0:
                    data.append({
                        'species': row.get('species', 'Unknown'),
                        'genus': row.get('genus', 'Unknown'),
                        'family': row.get('family', 'Unknown'),
                        'class': row.get('class', 'Unknown'),
                        'phylum': row.get('phylum', 'Unknown'),
                        'superkingdom': row.get('superkingdom', 'Unknown'),
                        'abundance': abundance,
                        'counts': float(row.get('estimated counts', 0))
                    })
    except Exception as e:
        print(f"Warning: Could not parse {file_path}: {e}", file=sys.stderr)
    return data

def parse_combined_abundance(file_path):
    data = {"samples": [], "taxa": []}
    try:
        with open(file_path, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            # Extract sample names (columns after the taxonomic ranks)
            tax_cols = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
            data["samples"] = [col for col in reader.fieldnames if col not in tax_cols]
            
            for row in reader:
                # Convert abundances to float, handle empty/missing as 0
                abundances = []
                for s in data["samples"]:
                    try:
                        val = float(row.get(s, 0) or 0)
                        abundances.append(val)
                    except ValueError:
                        abundances.append(0.0)
                
                if any(a > 0 for a in abundances):
                    taxon_entry = {col: row.get(col, 'Unknown') for col in tax_cols}
                    taxon_entry["abundances"] = abundances
                    data["taxa"].append(taxon_entry)
    except Exception as e:
        print(f"Warning: Could not parse combined results {file_path}: {e}", file=sys.stderr)
    return data

def main():
    if len(sys.argv) < 3:
        print("Usage: make_html_report.py <summary_tsv> <combined_tsv> [rel_abundance_files...]")
        sys.exit(1)

    summary_tsv = sys.argv[1]
    combined_tsv = sys.argv[2]
    rel_abundance_files = sys.argv[3:]

    # Parse summary counts
    summary_data = []
    headers = []
    try:
        with open(summary_tsv, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
            headers = reader.fieldnames
            for row in reader:
                summary_data.append(row)
    except Exception as e:
        print(f"Error: Could not parse {summary_tsv}: {e}", file=sys.stderr)

    # Parse combined results for heatmap
    heatmap_data = parse_combined_abundance(combined_tsv)
    
    # Natural sort heatmap samples and reorder abundances
    if heatmap_data["samples"]:
        original_samples = list(heatmap_data["samples"])
        sorted_samples = sorted(original_samples, key=natural_sort_key)
        idx_map = [original_samples.index(s) for s in sorted_samples]
        
        heatmap_data["samples"] = sorted_samples
        for taxon in heatmap_data["taxa"]:
            taxon["abundances"] = [taxon["abundances"][i] for i in idx_map]

    # Parse all rel-abundance files for bar plots
    all_samples = []
    for fpath in rel_abundance_files:
        sample_name = os.path.basename(fpath).replace('_rel-abundance.tsv', '').replace('.fastq_rel-abundance.tsv', '')
        data = parse_rel_abundance(fpath)
        if data:
            all_samples.append({'name': sample_name, 'data': data})
    
    # Natural sort all_samples by name
    all_samples.sort(key=lambda x: natural_sort_key(x['name']))

    table_header_html = "".join([f'<th class="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider sortable" data-key="{h.replace(" ", "_")}" onclick="sortSummaryTable({i})">{h}</th>' for i, h in enumerate(headers)]) if headers else ""
    
    table_rows_html = ""
    # Natural sort summary data by sample name
    summary_data.sort(key=lambda x: natural_sort_key(x.get('sample', '')))
    
    for row in summary_data:
        # Build data attributes for sorting using the raw keys
        data_attrs = " ".join([f'data-{h.replace(" ", "_")}="{row.get(h, "")}"' for h in headers])
        table_rows_html += f"<tr {data_attrs} class='hover:bg-gray-50'>"
        
        # Get basis for percentage text calculation
        total_reads = float(row.get('filtered_reads', row.get('raw_reads', 0)) or 1)
        if total_reads == 0: total_reads = 1
        
        for h in headers:
            val = row.get(h, "")
            display_val = val
            
            # Formatting logic matching nxf-alignment summary stats
            if h in ['mapped', 'unmapped', 'mapped_filtered', 'mapped_unclassified']:
                try:
                    count_val = float(val)
                    perc = (count_val / total_reads) * 100
                    display_val = f"{int(count_val):,} ({perc:.1f}%)"
                except: pass
            elif h in ['raw_reads', 'filtered_reads']:
                try: display_val = f"{int(float(val)):,}"
                except: pass
            elif h != 'sample':
                try: 
                    fval = float(val)
                    display_val = f"{int(fval):,}" if fval.is_integer() else f"{fval:.2f}"
                except: pass
            
            alignment = "text-right" if h != 'sample' else "text-left sample-col"
            # Consistent padding (10px approx px-4 py-2.5)
            table_rows_html += f'<td class="px-4 py-2.5 whitespace-nowrap text-sm text-gray-900 {alignment}">{display_val}</td>'
        table_rows_html += "</tr>"

    # Load external JavaScript from the same directory as the script
    js_template_path = os.path.join(os.path.dirname(__file__), 'report.js')
    try:
        with open(js_template_path, 'r') as f:
            js_content = f.read()
        js_content = js_content.replace('__SAMPLES_DATA__', json.dumps(all_samples))
        js_content = js_content.replace('__HEATMAP_DATA__', json.dumps(heatmap_data))
    except Exception as e:
        print(f"Warning: Could not load {js_template_path}: {e}", file=sys.stderr)
        js_content = f"// Error loading report.js: {e}"

    html_template = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Emu Abundance Report</title>
    <script src="https://cdn.tailwindcss.com"></script>
    <script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.2/dist/chart.umd.min.js"></script>
    <style>
        body { font-family: 'Inter', sans-serif; }
        .heatmap-cell { transition: all 0.2s; }
        .heatmap-cell:hover { transform: scale(1.1); z-index: 10; }
        
        /* Table Styles from nxf-alignment */
        #summaryTable { table-layout: fixed; width: 100%; }
        table { border-collapse: collapse; width: auto; min-width: 100%; }
        th, td { font-family: 'Courier New', monospace; border-bottom: 1px solid #e5e7eb; }
        th { font-family: 'Inter', sans-serif; }
        .sample-col { 
            font-weight: 600; 
            background: #f8fafc !important; 
            width: 200px; 
            min-width: 200px; 
            white-space: nowrap; 
            overflow: hidden; 
            text-overflow: ellipsis; 
        }
        
        /* Table Sorting Styles */
        th.sortable { cursor: pointer; user-select: none; position: relative; padding-right: 20px; }
        th.sortable:hover { background: #f3f4f6; }
        th.sortable::after { content: '⇅'; position: absolute; right: 5px; opacity: 0.3; font-size: 0.8em; }
        th.sortable.asc::after { content: '▲'; opacity: 1; }
        th.sortable.desc::after { content: '▼'; opacity: 1; }

        /* Collapsible Sections */
        .collapsible-section { margin-bottom: 2rem; border: 1px solid #e5e7eb; border-radius: 1rem; background: #ffffff; box-shadow: 0 4px 6px -1px rgb(0 0 0 / 0.1); }
        .collapsible-section summary { list-style: none; cursor: pointer; padding: 1.5rem 2rem; background: #ffffff; transition: background 0.2s; border-bottom: 1px solid #f3f4f6; }
        .collapsible-section summary:hover { background: #f9fafb; }
        .collapsible-section summary::-webkit-details-marker { display: none; }
        .collapsible-section summary h1 { display: flex; align-items: center; margin-bottom: 0; border-bottom: none; padding-bottom: 0; }
        .collapsible-section summary h1::before { 
            content: '▶'; display: inline-block; margin-right: 1rem; font-size: 0.7em; 
            transition: transform 0.2s; opacity: 0.3; 
        }
        .collapsible-section[open] summary h1::before { transform: rotate(90deg); }
        .collapsible-section .section-content { padding: 2rem; }
    </style>
</head>
<body class="bg-gray-100 min-h-screen p-6">
    <div class="max-w-7xl mx-auto space-y-8">
        <!-- Pipeline Summary -->
        <details class="collapsible-section" open>
            <summary><h1 class="text-xl font-bold text-gray-900">Pipeline Summary</h1></summary>
            <div class="section-content">
                <div class="overflow-x-auto border rounded-xl">
                    <table id="summaryTable" class="min-w-full divide-y divide-gray-200">
                        <thead class="bg-gray-50">
                            <tr>{{TABLE_HEADER}}</tr>
                        </thead>
                        <tbody class="bg-white divide-y divide-gray-200">{{TABLE_ROWS}}</tbody>
                    </table>
                </div>
            </div>
        </details>

        <!-- Global Controls -->
        <details class="collapsible-section" open>
            <summary><h1 class="text-xl font-bold text-gray-900">Report Controls</h1></summary>
            <div class="section-content">
                <div class="grid grid-cols-1 md:grid-cols-3 gap-4 bg-gray-50 p-4 rounded-xl">
                    <div>
                        <label class="block text-sm font-medium text-gray-700 mb-1">Group By</label>
                        <select id="rankSelect" class="w-full border-gray-300 rounded-md shadow-sm p-2 outline-none focus:ring-2 focus:ring-indigo-500">
                            <option value="superkingdom">Superkingdom</option>
                            <option value="phylum">Phylum</option>
                            <option value="class">Class</option>
                            <option value="family">Family</option>
                            <option value="genus" selected>Genus</option>
                            <option value="species">Species</option>
                        </select>
                    </div>
                    <div>
                        <label class="block text-sm font-medium text-gray-700 mb-1">Top N Taxa</label>
                        <input type="number" id="topNInput" value="10" min="1" max="100" class="w-full border-gray-300 rounded-md shadow-sm p-2 outline-none focus:ring-2 focus:ring-indigo-500">
                    </div>
                    <div>
                        <label class="block text-sm font-medium text-gray-700 mb-1">Select Samples</label>
                        <div class="relative">
                            <button id="dropdownButton" class="w-full flex justify-between items-center border border-gray-300 rounded-md shadow-sm p-2 bg-white text-sm outline-none focus:ring-2 focus:ring-indigo-500">
                                <span id="dropdownLabel">All Samples Selected</span>
                                <svg class="w-4 h-4 ml-2" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M19 9l-7 7-7-7"></path></svg>
                            </button>
                            <div id="sampleDropdown" class="absolute z-50 mt-1 w-full bg-white border border-gray-300 rounded-md shadow-lg hidden max-h-64 overflow-y-auto">
                                <div class="p-2 border-b border-gray-100 flex justify-between">
                                    <button onclick="toggleAllSamples(true)" class="text-xs text-indigo-600 hover:text-indigo-800 font-semibold">Select All</button>
                                    <button onclick="toggleAllSamples(false)" class="text-xs text-indigo-600 hover:text-indigo-800 font-semibold">Clear</button>
                                </div>
                                <div id="sampleSelector" class="p-1"></div>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </details>

        <!-- Abundance Plots -->
        <details class="collapsible-section" open>
            <summary><h1 class="text-xl font-bold text-gray-900">Taxonomic Distribution</h1></summary>
            <div class="section-content">
                <div id="chartContainer" class="relative min-h-[300px]">
                    <div id="chart-tooltip" class="absolute pointer-events-none z-50 p-3 bg-gray-800 text-white rounded shadow-xl opacity-0 transition-opacity duration-200 w-72 overflow-hidden"></div>
                    <canvas id="chartCanvas"></canvas>
                </div>
            </div>
        </details>

        <!-- Heatmap Section -->
        <details class="collapsible-section" open>
            <summary><h1 class="text-xl font-bold text-gray-900">Abundance Heatmap</h1></summary>
            <div class="section-content">
                <div class="overflow-x-auto border rounded-xl bg-gray-50">
                    <div id="heatmapContainer" class="p-4 inline-block min-w-full">
                        <!-- Heatmap will be generated here -->
                    </div>
                </div>
            </div>
        </details>
    </div>

    <script>
{{JS_CONTENT}}
    </script>
</body>
</html>
"""
    html_content = html_template.replace('{{TABLE_HEADER}}', table_header_html)
    html_content = html_content.replace('{{TABLE_ROWS}}', table_rows_html)
    html_content = html_content.replace('{{JS_CONTENT}}', js_content)
    
    with open('nxf-emu-report.html', 'w') as f:
        f.write(html_content)

if __name__ == "__main__":
    main()
