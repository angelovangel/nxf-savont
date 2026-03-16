#!/usr/bin/env python3
import sys
import os
import json
import csv
import re
import argparse
from datetime import datetime

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
                    try:
                        counts = float(row.get('estimated counts', 0))
                    except (ValueError, TypeError):
                        counts = 0
                    data.append({
                        'tax_id': row.get('tax_id', 'Unknown'),
                        'species': row.get('species', 'Unknown'),
                        'genus': row.get('genus', 'Unknown'),
                        'family': row.get('family', 'Unknown'),
                        'order': row.get('order', 'Unknown'),
                        'class': row.get('class', 'Unknown'),
                        'phylum': row.get('phylum', 'Unknown'),
                        'clade': row.get('clade', 'Unknown'),
                        'superkingdom': row.get('superkingdom', 'Unknown'),
                        'abundance': abundance,
                        'counts': counts
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
            tax_cols = ['tax_id', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'clade', 'superkingdom']
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

def render_wfinfo_block(title, info_list, add_top_border=False):
    if not info_list:
        return ""
    
    border_style = "border-top: 1px solid rgba(255, 255, 255, 0.2);" if add_top_border else ""
    html = f"""
    <details style="padding-top: 5px; margin-top: 5px; {border_style} text-align: left;">
      <summary style="font-size: 0.8em; cursor: pointer; color: white;">{title}</summary>
      <div style="padding-top: 5px; text-align: left; width: 100%;">
        <table style="width: 100%; border-collapse: collapse; margin-top: 5px; color: white; border: none; background: inherit;">
          <tbody>
    """
    
    # Process pairs since savont uses list of tuples, unlike alignment's list of dicts
    for i, (key, value) in enumerate(info_list):
        if i > 0:
            html += '<tr><td colspan="2" style="border-bottom: 1px solid rgba(255, 255, 255, 0.2); padding: 0;"></td></tr>'
        
        display_key = key.replace('_', ' ').title()
        html += f"""
        <tr>
          <td style="padding: 3px 10px 3px 0; font-weight: 500; border: none; background: inherit; font-family: 'Courier New', monospace; color: white; white-space: nowrap; font-size: 0.8em; width: 180px;">{display_key}:</td>
          <td style="padding: 3px 0; border: none; background: inherit; font-family: 'Courier New', monospace; color: white; font-size: 0.8em;">{value}</td>
        </tr>
        """
    
    html += """
          </tbody>
        </table>
      </div>
    </details>
    """
    return html

    return html

def main():
    parser = argparse.ArgumentParser(description='Generate Savont Abundance Report')
    parser.add_argument('--summary', required=True, help='Summary TSV file')
    parser.add_argument('--combined', required=True, help='Combined TSV file')
    parser.add_argument('--abundances', nargs='*', help='Relative abundance TSV files', default=[])
    parser.add_argument('--wfinfo', help='Optional CSV file with workflow properties', default=None)
    
    args = parser.parse_args()

    summary_tsv = args.summary
    combined_tsv = args.combined
    rel_abundance_files = args.abundances
    wfinfo_csv = args.wfinfo

    # Parse workflow info
    wf_info = []
    if wfinfo_csv and os.path.exists(wfinfo_csv):
        try:
            with open(wfinfo_csv, 'r') as f:
                reader = csv.reader(f)
                for row in reader:
                    if len(row) >= 2:
                        wf_info.append((row[0], row[1]))
        except Exception as e:
            print(f"Warning: Could not parse workflow info {wfinfo_csv}: {e}", file=sys.stderr)

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

    table_header_html = ""
    for i, h in enumerate(headers):
        align = "text-left" if h == 'sample' else "text-right pr-6"
        table_header_html += f'<th class="px-4 py-3 {align} text-xs font-medium text-gray-500 uppercase tracking-wider sortable" data-key="{h.replace(" ", "_")}" onclick="sortSummaryTable({i})">{h}</th>'
    
    table_rows_html = ""
    # Natural sort summary data by sample name
    summary_data.sort(key=lambda x: natural_sort_key(x.get('sample', '')))
    
    for row in summary_data:
        # Build data attributes for sorting using the raw keys
        data_attrs = " ".join([f'data-{h.replace(" ", "_")}="{row.get(h, "")}"' for h in headers])
        table_rows_html += f"<tr {data_attrs} class='hover:bg-gray-50'>"
        
        # Get basis for percentage text calculation
        total_reads = float(row.get('subsampled_reads', row.get('filtered_reads', row.get('raw_reads', 0))) or 1)
        if total_reads == 0: total_reads = 1
        
        for h in headers:
            val = row.get(h, "")
            display_val = val
            
            # Formatting logic matching nxf-alignment summary stats
            if h in ['raw_reads', 'filtered_reads', 'raw_n50']:
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
        js_content = js_content.replace('__SUMMARY_DATA__', json.dumps(summary_data))
    except Exception as e:
        print(f"Warning: Could not load {js_template_path}: {e}", file=sys.stderr)
        js_content = f"// Error loading report.js: {e}"

    # Render Workflow Info
    wf_info_block = render_wfinfo_block("Workflow Details", wf_info, add_top_border=True)

    html_template = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Savont Abundance Report</title>
    <script src="https://cdn.tailwindcss.com"></script>
    <script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.2/dist/chart.umd.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/chroma-js/2.4.2/chroma.min.js"></script>
    <script src="https://d3js.org/d3.v7.min.js"></script>
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

        /* Report Header */
        .report-header {
            background: #374151;
            color: white;
            padding: 8px 10px;
            margin-bottom: 2rem;
            border-radius: 1rem;
            box-shadow: 0 4px 6px -1px rgb(0 0 0 / 0.15);
            text-align: center;
        }
        .report-header-main {
            display: flex;
            justify-content: center;
            align-items: center;
            gap: 15px;
            margin-bottom: 2px;
        }
        .report-header h2 {
            font-size: 1.2em;
            margin-bottom: 0;
            font-weight: 700;
        }
        .report-datetime {
            opacity: 0.9;
            font-size: 0.8em;
            margin-bottom: 0;
            text-align: right;
        }
        .repo-link {
            color: white;
            opacity: 0.7;
            transition: opacity 0.2s;
            display: flex;
            align-items: center;
            text-decoration: none;
        }
        .repo-link:hover { opacity: 1; }
    </style>
</head>
<body class="bg-gray-100 min-h-screen p-6">
    <div class="max-w-7xl mx-auto space-y-8">
        <!-- Report Header -->
        <div class="report-header">
            <div class="report-header-main">
                <h2>NXF-SAVONT Report</h2>
                <a href="https://github.com/angelovangel/nxf-savont" class="repo-link" target="_blank" title="View Source on GitHub">
                    <svg height="22" viewBox="0 0 16 16" version="1.1" width="22" aria-hidden="true" fill="currentColor"><path d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z"></path></svg>
                </a>
            </div>
            <div class="report-datetime">{REPORT_DATETIME}</div>
            {{WF_INFO}}
        </div>
        <!-- Summary -->
        <details class="collapsible-section" open>
            <summary>
                <div class="flex justify-between items-center w-full pr-8">
                    <h1 class="text-xl text-gray-900">Summary</h1>
                    <button onclick="downloadSummaryCSV(event)" class="flex items-center gap-2 px-3 py-1.5 bg-indigo-50 text-indigo-700 rounded-lg hover:bg-indigo-100 transition-colors text-sm font-semibold">
                        <svg class="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M4 16v1a3 3 0 003 3h10a3 3 0 003-3v-1m-4-4l-4 4m0 0l-4-4m4 4V4"></path></svg>
                        Download CSV
                    </button>
                </div>
            </summary>
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
            <summary><h1 class="text-xl text-gray-900">Report Controls</h1></summary>
            <div class="section-content">
                <div class="grid grid-cols-1 md:grid-cols-4 gap-4 bg-gray-50 p-4 rounded-xl">
                    <div>
                        <label class="block text-sm font-medium text-gray-700 mb-1">Color Scheme</label>
                        <select id="colorSchemeSelect" class="w-full border-gray-300 rounded-md shadow-sm p-2 outline-none focus:ring-2 focus:ring-indigo-500">
                            <option value="Set1">Set1</option>
                            <option value="Dark2">Dark2</option>
                            <option value="Set2">Set2</option>
                            <option value="Accent">Accent</option>
                            <option value="Paired">Paired</option>
                            <option value="Pastel1">Pastel1</option>
                            <option value="Pastel2">Pastel2</option>
                        </select>
                    </div>
                    <div>
                        <label class="block text-sm font-medium text-gray-700 mb-1">Group By</label>
                        <select id="rankSelect" class="w-full border-gray-300 rounded-md shadow-sm p-2 outline-none focus:ring-2 focus:ring-indigo-500">
                            <option value="superkingdom">Superkingdom</option>
                            <option value="phylum">Phylum</option>
                            <option value="clade">Clade</option>
                            <option value="class">Class</option>
                            <option value="order">Order</option>
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
            <summary><h1 class="text-xl text-gray-900">Taxonomic Distribution</h1></summary>
            <div class="section-content">
                <div id="chartContainer" class="relative min-h-[300px]">
                    <div id="chart-tooltip" class="absolute pointer-events-none z-50 p-3 bg-gray-800 text-white rounded shadow-xl opacity-0 transition-opacity duration-200 w-72 overflow-hidden"></div>
                    <canvas id="chartCanvas"></canvas>
                </div>
            </div>
        </details>

        <!-- Heatmap Section -->
        <details class="collapsible-section" open>
            <summary>
                <div class="flex justify-between items-center w-full pr-8">
                    <h1 class="text-xl text-gray-900">Abundance Heatmap</h1>
                    <button onclick="downloadHeatmapCSV(event)" class="flex items-center gap-2 px-3 py-1.5 bg-indigo-50 text-indigo-700 rounded-lg hover:bg-indigo-100 transition-colors text-sm font-semibold">
                        <svg class="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M4 16v1a3 3 0 003 3h10a3 3 0 003-3v-1m-4-4l-4 4m0 0l-4-4m4 4V4"></path></svg>
                        Download CSV
                    </button>
                </div>
            </summary>
            <div class="section-content">
                <div class="overflow-x-auto border rounded-xl bg-gray-50">
                    <div id="heatmapContainer" class="p-4 inline-block min-w-full">
                        <!-- Heatmap will be generated here -->
                    </div>
                </div>
            </div>
        </details>

        <!-- Tree Section -->
        <details class="collapsible-section" open>
            <summary><h1 class="text-xl text-gray-900">Phylogenetic Tree</h1></summary>
            <div class="section-content">
                <div class="mb-4 flex items-center gap-4 bg-gray-50 p-4 rounded-xl">
                    <div class="flex-1">
                        <label class="block text-sm font-medium text-gray-700 mb-1">Select Sample for Tree View</label>
                        <select id="treeSampleSelect" class="w-full border-gray-300 rounded-md shadow-sm p-2 outline-none focus:ring-2 focus:ring-indigo-500" onchange="renderStandardTree()">
                            <!-- Options populated by JS -->
                        </select>
                    </div>
                    <div class="w-48">
                        <label class="block text-sm font-medium text-gray-700 mb-1">Abundance Cutoff</label>
                        <select id="cutoffSelect" class="w-full border-gray-300 rounded-md shadow-sm p-2 outline-none focus:ring-2 focus:ring-indigo-500" onchange="renderStandardTree()">
                            <option value="0">Show all</option>
                            <option value="0.01">> 1%</option>
                            <option value="0.001">> 0.1%</option>
                            <option value="0.0001">> 0.01%</option>
                        </select>
                    </div>
                </div>
                <div class="bg-white rounded-xl border border-gray-200 overflow-hidden relative min-h-[600px]" id="treeContainer">
                    <svg id="treeSvg" class="w-full" style="min-height: 600px;"></svg>
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
    html_content = html_content.replace('{{WF_INFO}}', wf_info_block)
    html_content = html_content.replace('{{JS_CONTENT}}', js_content)
    html_content = html_content.replace('{REPORT_DATETIME}', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    
    with open('nxf-savont-report.html', 'w') as f:
        f.write(html_content)

if __name__ == "__main__":
    main()
