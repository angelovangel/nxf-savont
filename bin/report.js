const samples = __SAMPLES_DATA__;
const heatmapData = __HEATMAP_DATA__;
const summaryData = __SUMMARY_DATA__;
const lineageData = typeof __LINEAGE_DATA__ !== 'undefined' ? __LINEAGE_DATA__ : {};
const rarefactionData = typeof __RAREFACTION_DATA__ !== 'undefined' ? __RAREFACTION_DATA__ : null;
let chartInstance = null;

// --- Dropdown Management ---
document.getElementById('dropdownButton').addEventListener('click', () => {
    document.getElementById('sampleDropdown').classList.toggle('hidden');
});

document.addEventListener('click', (e) => {
    const dropdown = document.getElementById('sampleDropdown');
    const button = document.getElementById('dropdownButton');
    if (!dropdown.contains(e.target) && !button.contains(e.target)) {
        dropdown.classList.add('hidden');
    }
});

function toggleAllSamples(checked) {
    document.querySelectorAll('#sampleSelector input').forEach(cb => cb.checked = checked);
    updateChart();
    renderHeatmap();
    if (typeof renderRarefaction === 'function') renderRarefaction();
}

function populateSampleSelector() {
    const container = document.getElementById('sampleSelector');
    samples.forEach((s, i) => {
        const div = document.createElement('div');
        div.className = 'flex items-center px-2 py-1 hover:bg-gray-100 rounded cursor-pointer';
        div.onclick = (e) => {
            if (e.target.tagName !== 'INPUT') {
                const cb = div.querySelector('input');
                cb.checked = !cb.checked;
                updateChart();
                renderHeatmap();
                if (typeof renderRarefaction === 'function') renderRarefaction();
            }
        };
        div.innerHTML = `<input type="checkbox" id="sample-${i}" value="${s.name}" checked class="mr-2 cursor-pointer">
                         <label for="sample-${i}" class="text-sm truncate w-full cursor-pointer">${s.name}</label>`;
        div.querySelector('input').addEventListener('change', (e) => {
            e.stopPropagation();
            updateChart();
            renderHeatmap();
            if (typeof renderRarefaction === 'function') renderRarefaction();
        });
        container.appendChild(div);
    });
}

function getSelectedSampleNames() {
    return Array.from(document.querySelectorAll('#sampleSelector input:checked')).map(cb => cb.value);
}

// --- Heatmap Logic ---
function getTaxonName(d, rank) {
    if (d.tax_id && lineageData[d.tax_id]) {
        const tkLineage = lineageData[d.tax_id];
        let match = tkLineage.find(l => l.rank === rank);
        if (rank === 'superkingdom' && !match) {
            match = tkLineage.find(l => l.rank === 'domain');
        }
        if (match) return match.name;
    }
    return d[rank] || 'Unknown';
}

function getHeatmapColor(val) {
    if (val <= 0) return { bg: '#f9fafb', text: 'black' };
    const scaled = Math.pow(val, 0.4);
    const hue = 60 * (1 - scaled); // From 60 (Yellow) to 0 (Red)
    const saturation = 20 + (scaled * 80); // From 20% to 100%
    const lightness = 98 - (scaled * 63);   // From 98% (pale) to 35% (deep red)
    return {
        bg: `hsl(${hue}, ${saturation}%, ${lightness}%)`,
        text: lightness < 55 ? 'white' : 'black'
    };
}

function renderHeatmap() {
    const rank = document.getElementById('rankSelect').value;
    const topN = parseInt(document.getElementById('topNInput').value) || 10;
    const container = document.getElementById('heatmapContainer');
    const selectedNames = getSelectedSampleNames();

    if (selectedNames.length === 0) {
        container.innerHTML = '<p class="text-gray-500 italic p-4">No samples selected</p>';
        return;
    }

    const sampleIndices = selectedNames.map(name => {
        let idx = heatmapData.samples.indexOf(name);
        if (idx === -1) {
            // Try fuzzy match: strip .fastq, .filtered, .tsv etc.
            const clean = (s) => s.replace(/(\.filtered|\.fastq|\.tsv|_rel-abundance\.tsv)$/g, '');
            const cleanName = clean(name);
            idx = heatmapData.samples.findIndex(s => clean(s) === cleanName);
        }
        return idx;
    }).filter(idx => idx !== -1);

    // Aggregate by selected rank
    const aggregated = {};
    const taxonTaxidSets = {}; // taxonName -> Set(taxids)
    const taxonLineages = {}; // taxonName -> lineage object

    const ranks = ['superkingdom', 'phylum', 'clade', 'class', 'order', 'family', 'genus', 'species'];
    const rankIndex = ranks.indexOf(rank);

    heatmapData.taxa.forEach(t => {
        const taxonName = getTaxonName(t, rank);
        if (!aggregated[taxonName]) {
            aggregated[taxonName] = new Array(heatmapData.samples.length).fill(0);
            taxonTaxidSets[taxonName] = new Set();

            // Capture lineage up to current rank
            let lineage = [];
            const taxonKitLineage = lineageData[t.tax_id];

            if (taxonKitLineage) {
                const rankIdx = taxonKitLineage.findIndex(l => l.name === taxonName);
                if (rankIdx !== -1) {
                    lineage = taxonKitLineage.slice(0, rankIdx + 1);
                } else {
                    lineage = taxonKitLineage;
                }
            } else {
                for (let i = 0; i <= rankIndex; i++) {
                    const r = ranks[i];
                    if (t[r] && t[r] !== 'Unknown') {
                        lineage.push({ rank: r, name: t[r] });
                    }
                }
            }
            taxonLineages[taxonName] = lineage;
        }
        t.abundances.forEach((val, i) => {
            aggregated[taxonName][i] += val;
        });
        if (rank === 'species' && t.tax_id && t.tax_id !== 'Unknown' && t.tax_id !== '0') {
            taxonTaxidSets[taxonName].add(t.tax_id);
        }
    });

    const allTaxa = Object.keys(aggregated).map(name => {
        const ids = Array.from(taxonTaxidSets[name]);
        return {
            name: name,
            sum: sampleIndices.reduce((acc, idx) => acc + (aggregated[name][idx] || 0), 0),
            abundances: aggregated[name],
            taxid: ids.length === 1 ? ids[0] : null,
            lineage: taxonLineages[name]
        };
    })
        .filter(t => t.sum > 0)
        .sort((a, b) => b.sum - a.sum);

    const sortedTaxa = allTaxa.slice(0, topN);

    // Calculate "Other" row
    if (allTaxa.length > topN) {
        const otherAbundances = new Array(heatmapData.samples.length).fill(0);
        allTaxa.slice(topN).forEach(t => {
            t.abundances.forEach((val, i) => {
                otherAbundances[i] += val;
            });
        });
        sortedTaxa.push({
            name: 'Other',
            sum: otherAbundances.reduce((acc, val) => acc + val, 0),
            abundances: otherAbundances
        });
    }

    let html = '<table class="border-collapse text-xs" style="table-layout: auto; min-width: max-content;">';
    // Header
    html += '<thead><tr><th class="p-2 text-left bg-white border sticky left-0 z-20">Taxon</th>';
    selectedNames.forEach(name => {
        html += `<th class="p-2 text-center bg-gray-100 border whitespace-nowrap">${name}</th>`;
    });
    html += '</tr></thead><tbody>';

    // Rows
    sortedTaxa.forEach(t => {
        const isOther = t.name === 'Other';
        const taxidAttr = t.taxid ? `data-taxid="${t.taxid}"` : '';
        const lineageAttr = t.lineage ? `data-lineage='${JSON.stringify(t.lineage).replace(/'/g, "&apos;")}'` : '';
        const hoverClass = (t.taxid || t.lineage) ? 'taxon-hover cursor-help underline decoration-dotted decoration-indigo-300 underline-offset-4' : '';

        html += `<tr><td class="p-2 font-medium ${isOther ? 'bg-gray-50 italic' : 'bg-white'} border sticky left-0 z-10 whitespace-nowrap shadow-sm ${hoverClass}" ${taxidAttr} ${lineageAttr}>${t.name}</td>`;
        sampleIndices.forEach(idx => {
            const val = t.abundances[idx] || 0;
            const { bg, text } = getHeatmapColor(val);
            html += `<td class="p-2 text-center border heatmap-cell" style="background-color: ${bg}; color: ${text}" title="${t.name} in ${heatmapData.samples[idx]}: ${(val * 100).toFixed(4)}%">
                            ${val > 0.0001 ? (val * 100).toFixed(2) + '%' : (val > 0 ? '<0.01%' : '-')}
                         </td>`;
        });
        html += '</tr>';
    });
    html += '</tbody></table>';
    container.innerHTML = html;

    // Attach hover listeners
    container.querySelectorAll('.taxon-hover').forEach(el => {
        el.onmouseenter = (e) => showTaxonTooltip(e, el.dataset.taxid, el.innerText, el.dataset.lineage);
        el.onmouseleave = hideTaxonTooltip;
    });
}

let tooltipTimeout = null;

function showTaxonTooltip(e, taxid, name, lineageJson) {
    const tooltip = document.getElementById('heatmap-tooltip');
    if (!tooltip) return;

    if (tooltipTimeout) {
        clearTimeout(tooltipTimeout);
        tooltipTimeout = null;
    }

    let lineageHtml = '';
    if (lineageJson) {
        try {
            const lineage = JSON.parse(lineageJson);
            lineageHtml = '<div class="px-4 py-2 bg-white border-b border-gray-100">';
            lineage.forEach((item, idx) => {
                const isLast = idx === lineage.length - 1;
                const indent = idx * 12;
                lineageHtml += `
                    <div style="padding-left: ${indent}px; display: flex; align-items: center; gap: 8px;">
                        ${idx > 0 ? '<span class="text-gray-300">└─</span>' : ''}
                        <span class="text-[11px] ${isLast ? 'font-bold text-gray-900' : 'text-gray-500'}">${item.name}</span>
                    </div>
                `;
            });
            lineageHtml += '</div>';
        } catch (err) {
            console.error("Error parsing lineage JSON", err);
        }
    }

    const ncbiLink = taxid ? `
        <div class="p-4 pt-2 bg-white">
            <a href="https://www.ncbi.nlm.nih.gov/datasets/taxonomy/${taxid}/" target="_blank" class="flex items-center gap-2 text-indigo-600 hover:text-indigo-800 text-xs font-semibold">
                <span>View NCBI Taxonomy</span>
                <svg class="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M10 6H6a2 2 0 00-2 2v10a2 2 0 002 2h10a2 2 0 002-2v-4M14 4h6m0 0v6m0-6L10 14"></path></svg>
            </a>
        </div>
    ` : '';

    tooltip.innerHTML = `
        <div id="tooltip-header" class="border-none bg-gray-50">
            <div>
                <div class="text-[10px] text-gray-400 uppercase tracking-wider font-bold leading-tight">Taxonomic Hierarchy</div>
                <div class="text-sm font-semibold italic text-gray-900">${name}</div>
                ${taxid ? `<div class="text-[10px] text-gray-400 font-mono">TaxID: ${taxid}</div>` : ''}
            </div>
        </div>
        ${lineageHtml}
        ${ncbiLink}
    `;

    tooltip.style.width = '280px';
    tooltip.style.display = 'flex';

    // Position handling - place to the right of the mouse
    // We use a small horizontal offset (10px) to prevent flickering 
    // but keep it close enough to maintain the "safe zone" for the mouse.
    let x = e.clientX + 10;
    let y = e.clientY - 20; // Slightly above to align better with the text line

    const rect = tooltip.getBoundingClientRect();

    // Check horizontal viewport bounds
    if (x + rect.width > window.innerWidth) {
        x = e.clientX - rect.width - 20; // Flip to left if no space on right
    }

    // Check vertical viewport bounds
    if (y + rect.height > window.innerHeight) {
        y = window.innerHeight - rect.height - 20;
    }

    if (x < 10) x = 10;
    if (y < 10) y = 10;

    tooltip.style.left = x + 'px';
    tooltip.style.top = y + 'px';

    tooltip.onmouseenter = () => {
        if (tooltipTimeout) {
            clearTimeout(tooltipTimeout);
            tooltipTimeout = null;
        }
    };
    tooltip.onmouseleave = hideTaxonTooltip;
}

function hideTaxonTooltip() {
    const tooltip = document.getElementById('heatmap-tooltip');
    if (tooltip) {
        if (tooltipTimeout) clearTimeout(tooltipTimeout);
        tooltipTimeout = setTimeout(() => {
            if (!tooltip.matches(':hover')) {
                tooltip.style.display = 'none';
            }
        }, 300); // Increased timeout for better feel
    }
}

// --- Chart Logic ---
let colorCache = {};
function getColor(name, index, total) {
    if (name === 'Other') return '#9ca3af';
    if (colorCache[name]) return colorCache[name];

    const scheme = document.getElementById('colorSchemeSelect').value;
    const colors = chroma.brewer[scheme] || chroma.scale(scheme).mode('lch').colors(Math.max(total, 10));

    const color = colors[index % colors.length];
    colorCache[name] = color;
    return color;
}

function externalTooltipHandler(context) {
    const tooltipEl = document.getElementById('chart-tooltip');
    const { chart, tooltip } = context;
    if (tooltip.opacity === 0) { tooltipEl.style.opacity = 0; return; }

    if (tooltip.body) {
        const datasetIndex = tooltip.dataPoints[0].datasetIndex;
        const dataset = chart.data.datasets[datasetIndex];
        const labels = chart.data.labels;
        let html = `<div class="font-bold border-b border-gray-600 mb-2 pb-1 text-sm mr-2" style="color: ${dataset.backgroundColor}">${dataset.label}</div><ul class="text-xs space-y-1">`;
        for (let i = 0; i < labels.length; i++) {
            const isHovered = i === tooltip.dataPoints[0].dataIndex ? 'bg-gray-700 rounded' : '';
            html += `<li class="flex justify-between items-center p-1 ${isHovered}"><span class="truncate mr-4">${labels[i]}</span><span class="font-mono">${dataset.data[i].toFixed(2)}%</span></li>`;
        }
        tooltipEl.innerHTML = html + '</ul>';
    }
    tooltipEl.style.opacity = 1;
    tooltipEl.style.left = chart.canvas.offsetLeft + 10 + 'px';
    tooltipEl.style.top = Math.min(chart.canvas.offsetTop + 10, chart.height - tooltipEl.offsetHeight - 10) + 'px';
}

function updateChart() {
    const rank = document.getElementById('rankSelect').value;
    const topN = parseInt(document.getElementById('topNInput').value) || 10;
    const selectedNames = getSelectedSampleNames();
    const selectedSamples = samples.filter(s => selectedNames.includes(s.name));

    document.getElementById('dropdownLabel').innerText = selectedNames.length === samples.length ? 'All Samples Selected' : (selectedNames.length === 0 ? 'No Samples Selected' : `${selectedNames.length} Samples Selected`);

    const container = document.getElementById('chartContainer');
    if (selectedSamples.length === 0) {
        if (chartInstance) chartInstance.destroy();
        chartInstance = null;
        container.style.height = '300px';
        return;
    }

    container.style.height = Math.max(300, selectedSamples.length * 40 + 100) + 'px';

    const globalTotals = {};
    const presence = {};
    selectedSamples.forEach(s => {
        s.data.forEach(d => {
            const taxon = getTaxonName(d, rank);
            globalTotals[taxon] = (globalTotals[taxon] || 0) + d.abundance;
            if (!presence[taxon]) presence[taxon] = new Set();
            if (d.abundance > 0) presence[taxon].add(s.name);
        });
    });

    const commonTaxa = Object.keys(globalTotals).filter(t => presence[t].size === selectedSamples.length);
    const topTaxa = commonTaxa.sort((a, b) => globalTotals[b] - globalTotals[a]).slice(0, topN);

    const datasets = topTaxa.map((taxon, i) => {
        const taxonData = selectedSamples.map(s => {
            const match = s.data.filter(d => getTaxonName(d, rank) === taxon);
            return { abundance: match.reduce((sum, d) => sum + d.abundance, 0) * 100, counts: Math.round(match.reduce((sum, d) => sum + d.counts, 0)) };
        });
        return { label: taxon, data: taxonData.map(d => d.abundance), counts: taxonData.map(d => d.counts), backgroundColor: getColor(taxon, i, topTaxa.length) };
    });

    const otherData = selectedSamples.map(s => {
        const other = s.data.filter(d => !topTaxa.includes(getTaxonName(d, rank)));
        return { abundance: other.reduce((sum, d) => sum + d.abundance, 0) * 100, counts: Math.round(other.reduce((sum, d) => sum + d.counts, 0)) };
    });

    if (otherData.some(v => v.abundance > 0)) {
        datasets.push({ label: 'Other', data: otherData.map(d => d.abundance), counts: otherData.map(d => d.counts), backgroundColor: '#9ca3af' });
    }

    if (chartInstance) chartInstance.destroy();
    chartInstance = new Chart(document.getElementById('chartCanvas'), {
        type: 'bar',
        data: { labels: selectedSamples.map(s => s.name), datasets },
        options: {
            indexAxis: 'y', responsive: true, maintainAspectRatio: false,
            barPercentage: 0.8, categoryPercentage: 1.0,
            scales: { x: { stacked: true, max: 100, title: { display: true, text: 'Relative Abundance (%)' } }, y: { stacked: true } },
            plugins: { legend: { position: 'right' }, tooltip: { enabled: false, external: externalTooltipHandler } }
        }
    });
}

document.getElementById('rankSelect').addEventListener('change', () => {
    updateChart();
    renderHeatmap();
});
document.getElementById('topNInput').addEventListener('change', () => {
    updateChart();
    renderHeatmap();
});
document.getElementById('colorSchemeSelect').addEventListener('change', () => {
    colorCache = {};
    updateChart();
    if (typeof renderRarefaction === 'function') renderRarefaction();
});

// --- Download Logic ---
function downloadCSV(filename, csvContent) {
    const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
    const link = document.createElement("a");
    if (link.download !== undefined) {
        const url = URL.createObjectURL(blob);
        link.setAttribute("href", url);
        link.setAttribute("download", filename);
        link.style.visibility = 'hidden';
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
    }
}

function downloadSummaryCSV(event) {
    if (event) event.stopPropagation();
    if (!summaryData || summaryData.length === 0) return;

    const headers = Object.keys(summaryData[0]);
    let csv = headers.join(',') + '\n';

    summaryData.forEach(row => {
        csv += headers.map(h => {
            let val = row[h] || "0";
            // Ensure read counts and related metrics are integers
            if (['raw_reads', 'filtered_reads', 'raw_n50'].includes(h)) {
                val = Math.round(parseFloat(val) || 0);
            }
            if (typeof val === 'string' && (val.includes(',') || val.includes('"') || val.includes('\n'))) {
                val = '"' + val.replace(/"/g, '""') + '"';
            }
            return val;
        }).join(',') + '\n';
    });

    downloadCSV('pipeline_summary.csv', csv);
}

function downloadHeatmapCSV(event) {
    if (event) event.stopPropagation();
    const rank = document.getElementById('rankSelect').value;
    const topN = parseInt(document.getElementById('topNInput').value) || 10;
    const selectedNames = getSelectedSampleNames();

    if (selectedNames.length === 0) {
        alert("No samples selected for download.");
        return;
    }

    const sampleIndices = selectedNames.map(name => {
        let idx = heatmapData.samples.indexOf(name);
        if (idx === -1) {
            const clean = (s) => s.replace(/(\.filtered|\.fastq|\.tsv|_rel-abundance\.tsv)$/g, '');
            const cleanName = clean(name);
            idx = heatmapData.samples.findIndex(s => clean(s) === cleanName);
        }
        return idx;
    }).filter(idx => idx !== -1);

    // Aggregate by selected rank
    const aggregated = {};
    heatmapData.taxa.forEach(t => {
        const taxonName = t[rank] || 'Unknown';
        if (!aggregated[taxonName]) {
            aggregated[taxonName] = new Array(heatmapData.samples.length).fill(0);
        }
        t.abundances.forEach((val, i) => {
            aggregated[taxonName][i] += val;
        });
    });

    const sortedTaxa = Object.keys(aggregated).map(name => ({
        name: name,
        sum: sampleIndices.reduce((acc, idx) => acc + (aggregated[name][idx] || 0), 0),
        abundances: aggregated[name]
    }))
        .filter(t => t.sum > 0)
        .sort((a, b) => b.sum - a.sum)
        .slice(0, topN);

    // Build CSV
    let csv = 'Taxon,' + selectedNames.join(',') + '\n';
    sortedTaxa.forEach(t => {
        let line = `"${t.name.replace(/"/g, '""')}"`;
        sampleIndices.forEach(idx => {
            line += ',' + (t.abundances[idx] || 0);
        });
        csv += line + '\n';
    });

    downloadCSV(`abundance_heatmap_${rank}.csv`, csv);
}

// --- Summary Table Sorting ---
let summarySortDirection = {};

function sortSummaryTable(columnIndex) {
    const table = document.getElementById('summaryTable');
    if (!table) return;
    const tbody = table.querySelector('tbody');
    const rows = Array.from(tbody.querySelectorAll('tr'));

    summarySortDirection[columnIndex] = summarySortDirection[columnIndex] === 'asc' ? 'desc' : 'asc';
    const isAsc = summarySortDirection[columnIndex] === 'asc';

    const headers = table.querySelectorAll('th.sortable');
    headers.forEach(h => h.classList.remove('asc', 'desc'));
    if (headers[columnIndex]) headers[columnIndex].classList.add(isAsc ? 'asc' : 'desc');

    rows.sort((a, b) => {
        const headerKey = headers[columnIndex].getAttribute('data-key');
        let aVal = a.getAttribute('data-' + headerKey) || "";
        let bVal = b.getAttribute('data-' + headerKey) || "";

        if (headerKey === 'sample' || isNaN(parseFloat(aVal)) || isNaN(parseFloat(bVal))) {
            return isAsc ? aVal.localeCompare(bVal, undefined, { numeric: true }) : bVal.localeCompare(aVal, undefined, { numeric: true });
        }

        return isAsc ? parseFloat(aVal) - parseFloat(bVal) : parseFloat(bVal) - parseFloat(aVal);
    });

    rows.forEach(row => tbody.appendChild(row));
}

// --- Tree Logic (D3) ---
function populateTreeSampleSelector() {
    const select = document.getElementById('treeSampleSelect');
    if (!select) return;
    select.innerHTML = '';
    samples.forEach(s => {
        const opt = document.createElement('option');
        opt.value = s.name;
        opt.innerText = s.name;
        select.appendChild(opt);
    });
}

function renderStandardTree() {
    const container = document.getElementById('treeContainer');
    if (!container) return;

    // Clear previous view
    container.innerHTML = '<svg id="treeSvg" class="w-full h-full"></svg>';
    const svgEl = document.getElementById('treeSvg');

    const sampleName = document.getElementById('treeSampleSelect').value;
    const cutoff = parseFloat(document.getElementById('cutoffSelect').value) || 0;
    const sample = samples.find(s => s.name === sampleName);

    if (!sample) return;

    // 1. Build hierarchical data
    const buildHierarchy = () => {
        const root = { name: "Life", children: [], tax_id: "root", rank: "root" };
        const fixedRanks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'];

        sample.data.forEach(t => {
            let current = root;
            if (t.abundance < cutoff) return; // Filter by cutoff

            let lineageArray = fixedRanks.map(targetRank => {
                let name = "Unknown";
                if (t.tax_id && lineageData[t.tax_id]) {
                    const tk = lineageData[t.tax_id];
                    let match = tk.find(l => l.rank === targetRank);
                    if (targetRank === 'superkingdom' && !match) {
                        match = tk.find(l => l.rank === 'domain');
                    }
                    if (match) name = match.name;
                } else {
                    name = t[targetRank] || "Unknown";
                }
                return { rank: targetRank, name: name };
            });

            lineageArray.forEach((item, i) => {
                const name = item.name;
                const rank = item.rank;
                if (!name) return;

                let child = current.children.find(c => c.name === name);
                if (!child) {
                    child = {
                        name: name,
                        rank: rank,
                        children: [],
                        tax_id: (i === fixedRanks.length - 1) ? t.tax_id : null
                    };
                    current.children.push(child);
                }

                if (i === fixedRanks.length - 1) {
                    child.value = (child.value || 0) + t.abundance;
                }
                current = child;
            });
        });
        return root;
    };

    const data = buildHierarchy();

    try {
        if (typeof d3 === 'undefined') {
            throw new Error("D3 library not loaded.");
        }

        const width = container.clientWidth || 1000;
        const height = 600;
        const margin = { top: 40, right: 250, bottom: 20, left: 40 };

        const svg = d3.select(svgEl)
            .attr("viewBox", [0, 0, width, height])
            .append("g")
            .attr("transform", `translate(${margin.left},${margin.top})`);

        const root = d3.hierarchy(data);
        root.sum(d => d.value || 0);

        // To align ranks, we need a fixed level step
        const levelStep = (width - margin.left - margin.right) / 7;

        const treeLayout = d3.cluster().size([height - margin.top - margin.bottom, width - margin.left - margin.right]);

        const radiusScale = d3.scaleSqrt()
            .domain([0, 1])
            .range([3, 22]);

        const phylumColorScale = d3.scaleOrdinal(d3.schemeCategory10);
        const defaultColor = "#94a3b8";

        const getPhylumColor = (d) => {
            let curr = d;
            while (curr && curr.data.rank !== 'phylum') curr = curr.parent;
            return curr ? phylumColorScale(curr.data.name) : defaultColor;
        };

        // Recursive expansion helper
        function expandToLevel(node, targetDepth) {
            if (node.depth < targetDepth) {
                if (node._children) {
                    node.children = node._children;
                    node._children = null;
                }
            } else if (node.depth >= targetDepth) {
                if (node.children) {
                    node._children = node.children;
                    node.children = null;
                }
            }
            if (node.children) {
                node.children.forEach(c => expandToLevel(c, targetDepth));
            } else if (node._children) {
                node._children.forEach(c => expandToLevel(c, targetDepth));
            }
        }

        // Draw rank headers and guide lines
        const displayRanks = ['Superkingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'];
        displayRanks.forEach((rank, i) => {
            const x = (i + 1) * levelStep;
            const targetDepth = i + 1;

            // Vertical dotted line
            svg.append("line")
                .attr("x1", x).attr("x2", x)
                .attr("y1", -20).attr("y2", height - margin.top - margin.bottom)
                .attr("stroke", "#e2e8f0")
                .attr("stroke-width", 1)
                .attr("stroke-dasharray", "4,4");

            // Rank label
            svg.append("text")
                .attr("x", x)
                .attr("y", -10)
                .attr("text-anchor", "middle")
                .attr("class", "text-[10px] font-bold text-gray-400 uppercase tracking-widest cursor-pointer hover:text-indigo-600 transition-colors")
                .style("cursor", "pointer")
                .text(rank)
                .on("click", () => {
                    expandToLevel(root, targetDepth);
                    update(root);
                });
        });

        let nodeCount = 0;

        function update(source) {
            treeLayout(root);
            const nodes = root.descendants();
            const links = root.links();

            // Align by rank
            nodes.forEach(d => {
                d.y = d.depth * levelStep;
            });

            // --- Links ---
            const link = svg.selectAll(".link")
                .data(links, d => d.target.id || (d.target.id = ++nodeCount));

            const linkEnter = link.enter()
                .append("path")
                .attr("class", "link")
                .attr("fill", "none")
                .attr("stroke", d => getPhylumColor(d.target))
                .attr("stroke-opacity", 0.4)
                .attr("stroke-width", 1.5)
                .attr("d", d => {
                    const o = { x: source.x, y: source.y };
                    return d3.linkHorizontal().x(d => d.y).y(d => d.x)({ source: o, target: o });
                });

            const linkUpdate = linkEnter.merge(link);
            linkUpdate.transition().duration(500)
                .attr("d", d3.linkHorizontal().x(d => d.y).y(d => d.x))
                .attr("stroke", d => getPhylumColor(d.target));

            link.exit().transition().duration(500)
                .attr("d", d => {
                    const o = { x: source.x, y: source.y };
                    return d3.linkHorizontal().x(d => d.y).y(d => d.x)({ source: o, target: o });
                })
                .remove();

            // --- Nodes ---
            const node = svg.selectAll(".node")
                .data(nodes, d => d.id || (d.id = ++nodeCount));

            const nodeEnter = node.enter()
                .append("g")
                .attr("class", "node")
                .attr("transform", d => `translate(${source.y},${source.x})`)
                .style("cursor", "pointer")
                .on("click", (event, d) => {
                    if (d.children) {
                        d._children = d.children;
                        d.children = null;
                    } else {
                        d.children = d._children;
                        d._children = null;
                    }
                    update(d);
                });

            // Hitbox for easy clicking
            nodeEnter.append("circle")
                .attr("r", 20)
                .attr("fill", "transparent");

            // Visible dot
            nodeEnter.append("circle")
                .attr("class", "visible-dot")
                .attr("r", 0)
                .attr("fill", d => getPhylumColor(d))
                .attr("stroke", "#fff")
                .attr("stroke-width", 1.5);

            nodeEnter.append("text")
                .attr("dy", "0.31em")
                .style("font-size", "11px")
                .style("fill-opacity", 0)
                .each(function () {
                    const el = d3.select(this);
                    el.append("tspan").attr("class", "name-tspan");
                    el.append("tspan").attr("class", "pct-tspan").style("font-weight", "bold").style("font-size", "12px").style("fill", "#64748b");
                });

            const nodeUpdate = nodeEnter.merge(node);

            nodeUpdate.transition().duration(500)
                .attr("transform", d => `translate(${d.y},${d.x})`);

            nodeUpdate.select(".visible-dot")
                .transition().duration(500)
                .attr("r", d => radiusScale(d.value))
                .attr("fill", d => getPhylumColor(d))
                .style("fill-opacity", d => (d._children || !d.children) ? 1 : 0.6); // Slightly fade expanded parent nodes

            const textUpdate = nodeUpdate.select("text");

            textUpdate.transition().duration(500)
                .style("fill-opacity", 1)
                .attr("x", d => d.children || d._children ? -(radiusScale(d.value) + 8) : (radiusScale(d.value) + 8))
                .attr("text-anchor", d => d.children || d._children ? "end" : "start");

            textUpdate.select(".name-tspan")
                .attr("x", d => d.children || d._children ? -(radiusScale(d.value) + 8) : (radiusScale(d.value) + 8))
                .attr("dy", "-0.1em")
                .text(d => d.data.name);

            textUpdate.select(".pct-tspan")
                .attr("x", d => d.children || d._children ? -(radiusScale(d.value) + 8) : (radiusScale(d.value) + 8))
                .attr("dy", "1.3em")
                .text(d => {
                    const pct = d.value * 100;
                    return pct < 0.1 && pct > 0 ? `${pct.toFixed(2)}%` : `${pct.toFixed(1)}%`;
                });

            nodeUpdate.select("title").remove();
            nodeUpdate.append("title")
                .text(d => `${d.data.name}\nAbundance: ${(d.value * 100).toFixed(4)}%`);

            const nodeExit = node.exit().transition().duration(500)
                .attr("transform", d => `translate(${source.y},${source.x})`)
                .remove();

            nodeExit.select("circle").attr("r", 0);
            nodeExit.select("text").style("fill-opacity", 0);
        }

        // Initialize with a manageable expansion (≤ 15 nodes at the deepest possible level)
        let bestDepth = 0;
        for (let d = 1; d <= 7; d++) {
            const count = root.descendants().filter(n => n.depth === d).length;
            if (count > 15) break;
            bestDepth = d;
        }

        if (bestDepth > 0) {
            expandToLevel(root, bestDepth);
        } else {
            // If even Superkingdom has > 10, just show Superkingdom collapsed
            expandToLevel(root, 1);
        }

        update(root);

    } catch (e) {
        console.error("Tree rendering error:", e);
        container.innerHTML += `<p class="text-red-500 p-4">Error rendering tree: ${e.message}</p>`;
    }
}

// --- Rarefaction Curve Logic ---
function renderRarefaction() {
    if (!rarefactionData || Object.keys(rarefactionData).length === 0) return;
    
    const RAW = rarefactionData;
    const samplesRaw = Object.keys(RAW);
    const container = document.getElementById('rarefactionChartContainer');
    const svgEl = document.getElementById('rarefactionSvg');
    if (!svgEl) return;
    
    const selectedNames = getSelectedSampleNames();
    const activeSamples = samplesRaw.filter(s => selectedNames.includes(s));
    
    if (activeSamples.length === 0) {
        svgEl.innerHTML = '';
        const msg = document.createElementNS("http://www.w3.org/2000/svg", "text");
        msg.setAttribute("x", "20");
        msg.setAttribute("y", "30");
        msg.setAttribute("fill", "#666");
        msg.textContent = "No samples selected.";
        svgEl.appendChild(msg);
        return;
    }
    
    svgEl.innerHTML = '';
    
    const margin = {top:20,right:150,bottom:50,left:65};
    const W = container.clientWidth || 820;
    const H = 440;
    const iW = W - margin.left - margin.right;
    const iH = H - margin.top  - margin.bottom;
    
    const svg = d3.select(svgEl).attr("viewBox", [0, 0, W, H]);
    const g   = svg.append("g").attr("transform",`translate(${margin.left},${margin.top})`);
    
    const xMax = d3.max(activeSamples.flatMap(s=>RAW[s].map(d=>d.depth)));
    const yMax = d3.max(activeSamples.flatMap(s=>RAW[s].map(d=>d.mean+d.std)));
    
    const x = d3.scaleLinear().domain([0,xMax*1.02]).range([0,iW]);
    const y = d3.scaleLinear().domain([0,yMax*1.1]).range([iH,0]);
    
    g.append("g").attr("class","grid").attr("transform",`translate(0,${iH})`)
      .call(d3.axisBottom(x).tickSize(-iH).tickFormat("").tickSizeOuter(0))
      .selectAll("line").attr("stroke", "#e8e8e8").attr("stroke-dasharray", "3,3");
      
    g.append("g").attr("class","grid")
      .call(d3.axisLeft(y).tickSize(-iW).tickFormat("").tickSizeOuter(0))
      .selectAll("line").attr("stroke", "#e8e8e8").attr("stroke-dasharray", "3,3");
      
    g.selectAll(".domain").remove();
    
    g.append("g").attr("class","axis").attr("transform",`translate(0,${iH})`)
      .call(d3.axisBottom(x).tickFormat(d3.format(",d")))
      .append("text").attr("x",iW/2).attr("y",40)
        .attr("fill","#333").attr("text-anchor","middle").style("font-size","12px")
        .text("Sequencing Depth (reads)");
        
    g.append("g").attr("class","axis").call(d3.axisLeft(y).ticks(6))
      .append("text").attr("transform","rotate(-90)").attr("x",-iH/2).attr("y",-50)
        .attr("fill","#333").attr("text-anchor","middle").style("font-size","12px")
        .text("Observed Species (Richness)");
        
    const areaGen = d3.area().x(d=>x(d.depth)).y0(d=>y(Math.max(0,d.mean-d.std))).y1(d=>y(d.mean+d.std));
    const lineGen = d3.line().x(d=>x(d.depth)).y(d=>y(d.mean));
    
    const sampleGroups = {};
    const muted = {};
    const tb = document.getElementById('toggleBand');
    const showBand = tb ? tb.checked : true;
    
    const tooltip = document.getElementById('rarefactionTooltip');
    
    activeSamples.forEach((s, i) => {
      const grp = g.append("g");
      sampleGroups[s] = grp;
      muted[s] = false;
      const c = getColor("sample_" + s, i, activeSamples.length);
      
      grp.append("path").datum(RAW[s]).attr("class","band").attr("fill",c)
         .style("opacity", 0.18).style("display", showBand ? null : "none")
         .attr("d",areaGen);
      grp.append("path").datum(RAW[s]).attr("class","line").attr("stroke",c)
         .attr("fill", "none").attr("stroke-width", 2.5)
         .attr("d",lineGen);
         
      grp.selectAll(".dot").data(RAW[s]).enter().append("circle").attr("class","dot")
        .attr("cx",d=>x(d.depth)).attr("cy",d=>y(d.mean)).attr("fill",c)
        .attr("r", 5).style("opacity", 0).style("cursor", "pointer")
        .on("mousemove",(event,d)=>{
          d3.select(event.currentTarget).style("opacity", 1);
          tooltip.style.display="block";
          tooltip.style.opacity = 1;
          const rect = container.getBoundingClientRect();
          let tx = event.clientX - rect.left + 14;
          let ty = event.clientY - rect.top - 14;
          tooltip.style.left= tx + "px";
          tooltip.style.top= ty + "px";
          tooltip.innerHTML=`<strong>${s}</strong><br>Depth: <b>${d3.format(",")(d.depth)}</b><br>Richness: <b>${d.mean.toFixed(2)}</b><br>&plusmn;SD: ${d.std.toFixed(2)}`;
        })
        .on("mouseleave",(event)=>{ 
          d3.select(event.currentTarget).style("opacity", 0);
          tooltip.style.opacity = 0; 
          tooltip.style.display="none"; 
        });
    });
    
    const legendG = svg.append("g").attr("transform",`translate(${margin.left+iW+18},${margin.top+10})`);
    legendG.append("text").text("Sample")
      .attr("font-size","12px").attr("font-weight","600").attr("fill","#333").attr("dy","0");
      
    activeSamples.forEach((s, i) => {
      const c = getColor("sample_" + s, i, activeSamples.length);
      const item = legendG.append("g").attr("class","leg-item")
        .attr("transform",`translate(0,${20 + i*22})`)
        .style("cursor","pointer")
        .on("click",()=>{
          muted[s]=!muted[s];
          item.style("opacity",muted[s]?0.3:1);
          sampleGroups[s].style("opacity",muted[s]?0:1).style("pointer-events",muted[s]?"none":"all");
        });
      item.append("line").attr("x1",0).attr("x2",22).attr("y1",0).attr("y2",0)
        .attr("stroke",c).attr("stroke-width",2.5);
      item.append("text").text(s)
        .attr("x",28).attr("dy","0.35em").attr("font-size","11px").attr("fill","#444");
    });
    
    if (tb) {
        tb.onchange = function() { 
            svg.selectAll(".band").style("display", this.checked ? null : "none"); 
        };
    }
}

// Initialize
window.onload = () => {
    populateSampleSelector();
    populateTreeSampleSelector();
    updateChart();
    renderHeatmap();
    renderStandardTree();
    if (typeof renderRarefaction === 'function') renderRarefaction();

    // Initial sort for summary table
    setTimeout(() => sortSummaryTable(0), 100);
};
