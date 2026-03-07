#!/usr/bin/env python3
import sys
import os
import json
import csv

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

def main():
    if len(sys.argv) < 2:
        print("Usage: make_html_report.py <summary_tsv> [rel_abundance_files...]")
        sys.exit(1)

    summary_tsv = sys.argv[1]
    rel_abundance_files = sys.argv[2:]

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

    # Parse all rel-abundance files
    all_samples = []
    for fpath in rel_abundance_files:
        sample_name = os.path.basename(fpath).replace('_rel-abundance.tsv', '').replace('.fastq_rel-abundance.tsv', '')
        data = parse_rel_abundance(fpath)
        if data:
            all_samples.append({'name': sample_name, 'data': data})

    table_header_html = "".join([f'<th class="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">{h}</th>' for h in headers]) if headers else ""
    
    table_rows_html = ""
    for row in summary_data:
        table_rows_html += "<tr>"
        
        # Get basis for percentage (filtered_reads or raw_reads)
        total_reads = float(row.get('filtered_reads', row.get('raw_reads', 0)) or 1)
        if total_reads == 0: total_reads = 1
        
        for h in headers:
            val = row.get(h, "")
            display_val = val
            
            # Format columns that are counts but not the basis columns
            if h in ['mapped', 'unmapped', 'mapped_filtered', 'mapped_unclassified']:
                try:
                    count_val = float(val)
                    perc = (count_val / total_reads) * 100
                    display_val = f"{int(count_val)} ({perc:.1f}%)"
                except:
                    pass
            elif h in ['raw_reads', 'filtered_reads']:
                try:
                    display_val = f"{int(float(val))}"
                except:
                    pass
            
            table_rows_html += f'<td class="px-6 py-4 whitespace-nowrap text-sm text-gray-500">{display_val}</td>'
        table_rows_html += "</tr>"

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Emu Abundance Report</title>
    <script src="https://cdn.tailwindcss.com"></script>
    <script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.2/dist/chart.umd.min.js"></script>
    <style>body {{ font-family: 'Inter', sans-serif; }}</style>
</head>
<body class="bg-gray-100 min-h-screen p-6">
    <div class="max-w-6xl mx-auto space-y-6">
        <div class="bg-white shadow-xl rounded-2xl p-8">
            <h1 class="text-3xl font-bold text-gray-900 mb-6 border-b pb-4">Pipeline Summary</h1>
            <div class="overflow-x-auto">
                <table class="min-w-full divide-y divide-gray-200">
                    <thead class="bg-gray-50">
                        <tr>{table_header_html}</tr>
                    </thead>
                    <tbody class="bg-white divide-y divide-gray-200">{table_rows_html}</tbody>
                </table>
            </div>
        </div>

        <div class="bg-white shadow-xl rounded-2xl p-8">
            <h1 class="text-3xl font-bold text-gray-900 mb-6 border-b pb-4">Microbial Abundance</h1>
            <div class="grid grid-cols-1 md:grid-cols-3 gap-4 mb-8 bg-gray-50 p-4 rounded-lg">
                <div>
                    <label class="block text-sm font-medium text-gray-700 mb-1">Group By</label>
                    <select id="rankSelect" class="w-full border-gray-300 rounded-md shadow-sm p-2">
                        <option value="superkingdom">Superkingdom</option>
                        <option value="phylum">Phylum</option>
                        <option value="class">Class</option>
                        <option value="family">Family</option>
                        <option value="genus" selected>Genus</option>
                        <option value="species">Species</option>
                    </select>
                </div>
                <div>
                    <label class="block text-sm font-medium text-gray-700 mb-1">Top N</label>
                    <input type="number" id="topNInput" value="10" min="1" max="30" class="w-full border-gray-300 rounded-md shadow-sm p-2">
                </div>
                <div>
                    <label class="block text-sm font-medium text-gray-700 mb-1">Select Samples</label>
                    <div class="relative">
                        <button id="dropdownButton" class="w-full flex justify-between items-center border border-gray-300 rounded-md shadow-sm p-2 bg-white text-sm">
                            <span id="dropdownLabel">All Samples Selected</span>
                            <svg class="w-4 h-4 ml-2" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M19 9l-7 7-7-7"></path></svg>
                        </button>
                        <div id="sampleDropdown" class="absolute z-50 mt-1 w-full bg-white border border-gray-300 rounded-md shadow-lg hidden max-h-64 overflow-y-auto">
                            <div class="p-2 border-b border-gray-100 flex justify-between">
                                <button onclick="toggleAllSamples(true)" class="text-xs text-indigo-600 hover:text-indigo-800 font-semibold">Select All</button>
                                <button onclick="toggleAllSamples(false)" class="text-xs text-indigo-600 hover:text-indigo-800 font-semibold">Clear</button>
                            </div>
                            <div id="sampleSelector" class="p-1">
                                <!-- Checkboxes will be inserted here -->
                            </div>
                        </div>
                    </div>
                </div>
            </div>

            <div id="chartContainer" class="relative min-h-[300px]">
                <div id="chart-tooltip" class="absolute pointer-events-none z-50 p-3 bg-gray-800 text-white rounded shadow-xl opacity-0 transition-opacity duration-200 w-72 overflow-hidden"></div>
                <canvas id="chartCanvas"></canvas>
            </div>
        </div>
    </div>

    <script>
        const samples = {json.dumps(all_samples)};
        let chartInstance = null;

        // Toggle dropdown
        document.getElementById('dropdownButton').addEventListener('click', () => {{
            document.getElementById('sampleDropdown').classList.toggle('hidden');
        }});

        // Close dropdown when clicking outside
        document.addEventListener('click', (e) => {{
            const dropdown = document.getElementById('sampleDropdown');
            const button = document.getElementById('dropdownButton');
            if (!dropdown.contains(e.target) && !button.contains(e.target)) {{
                dropdown.classList.add('hidden');
            }}
        }});

        function toggleAllSamples(checked) {{
            document.querySelectorAll('#sampleSelector input').forEach(cb => {{
                cb.checked = checked;
            }});
            updateChart();
        }}

        function populateSampleSelector() {{
            const container = document.getElementById('sampleSelector');
            samples.forEach((s, i) => {{
                const div = document.createElement('div');
                div.className = 'flex items-center px-2 py-1 hover:bg-gray-100 rounded cursor-pointer';
                div.onclick = (e) => {{
                    if (e.target.tagName !== 'INPUT') {{
                        const cb = div.querySelector('input');
                        cb.checked = !cb.checked;
                        updateChart();
                    }}
                }};
                div.innerHTML = `
                    <input type="checkbox" id="sample-${{i}}" value="${{s.name}}" checked class="mr-2 cursor-pointer">
                    <label for="sample-${{i}}" class="text-sm truncate w-full cursor-pointer">${{s.name}}</label>
                `;
                div.querySelector('input').addEventListener('change', (e) => {{
                    e.stopPropagation();
                    updateChart();
                }});
                container.appendChild(div);
            }});
        }}

        function getSelectedSamples() {{
            const selectedNames = Array.from(document.querySelectorAll('#sampleSelector input:checked')).map(cb => cb.value);
            return samples.filter(s => selectedNames.includes(s.name));
        }}

        function updateDropdownLabel(selectedCount) {{
            const label = document.getElementById('dropdownLabel');
            if (selectedCount === samples.length) {{
                label.innerText = 'All Samples Selected';
            }} else if (selectedCount === 0) {{
                label.innerText = 'No Samples Selected';
            }} else {{
                label.innerText = `${{selectedCount}} Samples Selected`;
            }}
        }}

        function getColor(name, index) {{
            if (name.startsWith('Other')) return '#9ca3af';
            const colors = ['#4f46e5', '#ef4444', '#10b981', '#f59e0b', '#3b82f6', '#8b5cf6', '#ec4899', '#6366f1', '#14b8a6', '#f97316'];
            return colors[index % colors.length];
        }}

        function externalTooltipHandler(context) {{
            const tooltipEl = document.getElementById('chart-tooltip');
            const {{chart, tooltip}} = context;

            if (tooltip.opacity === 0) {{
                tooltipEl.style.opacity = 0;
                return;
            }}

            if (tooltip.body) {{
                const datasetIndex = tooltip.dataPoints[0].datasetIndex;
                const dataset = chart.data.datasets[datasetIndex];
                const labels = chart.data.labels;
                
                let html = `<div class="font-bold border-b border-gray-600 mb-2 pb-1 text-sm mr-2" style="color: ${{dataset.backgroundColor}}">${{dataset.label}}</div>`;
                html += '<ul class="text-xs space-y-1">';
                
                for (let i = 0; i < labels.length; i++) {{
                    const abundance = dataset.data[i].toFixed(2);
                    const counts = dataset.counts[i].toLocaleString();
                    const isHovered = i === tooltip.dataPoints[0].dataIndex ? 'bg-gray-700 rounded' : '';
                    html += `<li class="flex justify-between items-center p-1 ${{isHovered}}">
                        <span class="truncate mr-4">${{labels[i]}}</span>
                        <span class="font-mono flex-shrink-0">${{abundance}}% (${{counts}})</span>
                    </li>`;
                }}
                html += '</ul>';
                tooltipEl.innerHTML = html;
            }}

            const {{offsetLeft: positionX, offsetTop: positionY}} = chart.canvas;
            tooltipEl.style.opacity = 1;
            tooltipEl.style.left = positionX + 10 + 'px';
            tooltipEl.style.top = Math.min(positionY + 10, chart.height - tooltipEl.offsetHeight - 10) + 'px';
        }}

        function updateChart() {{
            const rank = document.getElementById('rankSelect').value;
            const topN = parseInt(document.getElementById('topNInput').value) || 10;
            
            const selectedSamples = getSelectedSamples();
            updateDropdownLabel(selectedSamples.length);

            const container = document.getElementById('chartContainer');
            if (selectedSamples.length === 0) {{
                if (chartInstance) chartInstance.destroy();
                chartInstance = null;
                container.style.height = '300px';
                return;
            }}

            // Dynamically set height based on sample count to keep bars narrow but packed
            const targetHeight = Math.max(300, selectedSamples.length * 35 + 100);
            container.style.height = targetHeight + 'px';

            const globalTotals = {{}};
            const presence = {{}}; // taxon -> Set of sample names
            
            selectedSamples.forEach(s => {{
                s.data.forEach(d => {{
                    const taxon = d[rank] || 'Unknown';
                    globalTotals[taxon] = (globalTotals[taxon] || 0) + d.abundance;
                    if (!presence[taxon]) presence[taxon] = new Set();
                    if (d.abundance > 0) presence[taxon].add(s.name);
                }});
            }});

            // Filter for taxa present in EVERY selected sample
            const commonTaxa = Object.keys(globalTotals).filter(t => presence[t].size === selectedSamples.length);

            const topTaxa = commonTaxa
                .sort((a, b) => globalTotals[b] - globalTotals[a])
                .slice(0, topN);

            const datasets = topTaxa.map((taxon, i) => {{
                const taxonData = selectedSamples.map(s => {{
                    const match = s.data.filter(d => d[rank] === taxon);
                    return {{
                        abundance: match.reduce((sum, d) => sum + d.abundance, 0) * 100,
                        counts: Math.round(match.reduce((sum, d) => sum + d.counts, 0))
                    }};
                }});
                return {{
                    label: taxon,
                    data: taxonData.map(d => d.abundance),
                    counts: taxonData.map(d => d.counts),
                    backgroundColor: getColor(taxon, i)
                }};
            }});

            const otherData = selectedSamples.map(s => {{
                const other = s.data.filter(d => !topTaxa.includes(d[rank]));
                return {{
                    abundance: other.reduce((sum, d) => sum + d.abundance, 0) * 100,
                    counts: Math.round(other.reduce((sum, d) => sum + d.counts, 0))
                }};
            }});
            
            if (otherData.some(v => v.abundance > 0)) {{
                datasets.push({{
                    label: 'Other',
                    data: otherData.map(d => d.abundance),
                    counts: otherData.map(d => d.counts),
                    backgroundColor: '#9ca3af'
                }});
            }}

            if (chartInstance) chartInstance.destroy();
            chartInstance = new Chart(document.getElementById('chartCanvas'), {{
                type: 'bar',
                data: {{ labels: selectedSamples.map(s => s.name), datasets: datasets }},
                options: {{
                    indexAxis: 'y',
                    responsive: true,
                    maintainAspectRatio: false,
                    barPercentage: 0.8,
                    categoryPercentage: 1.0,
                    scales: {{ x: {{ stacked: true, max: 100, title: {{ display: true, text: 'Relative Abundance (%)' }} }}, y: {{ stacked: true }} }},
                    plugins: {{ 
                        legend: {{ position: 'right' }},
                        tooltip: {{ enabled: false, external: externalTooltipHandler }}
                    }}
                }}
            }});
        }}

        document.getElementById('rankSelect').addEventListener('change', updateChart);
        document.getElementById('topNInput').addEventListener('change', updateChart);
        populateSampleSelector();
        updateChart();
    </script>
</body>
</html>
"""
    with open('report.html', 'w') as f:
        f.write(html_content)

if __name__ == "__main__":
    main()
