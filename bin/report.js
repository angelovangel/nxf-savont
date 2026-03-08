const samples = __SAMPLES_DATA__;
const heatmapData = __HEATMAP_DATA__;
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
            }
        };
        div.innerHTML = `<input type="checkbox" id="sample-${i}" value="${s.name}" checked class="mr-2 cursor-pointer">
                         <label for="sample-${i}" class="text-sm truncate w-full cursor-pointer">${s.name}</label>`;
        div.querySelector('input').addEventListener('change', (e) => {
            e.stopPropagation();
            updateChart();
            renderHeatmap();
        });
        container.appendChild(div);
    });
}

function getSelectedSampleNames() {
    return Array.from(document.querySelectorAll('#sampleSelector input:checked')).map(cb => cb.value);
}

// --- Heatmap Logic ---
function getHeatmapColor(val) {
    if (val <= 0) return { bg: '#f9fafb', text: 'black' };
    const scaled = Math.pow(val, 0.4);
    const hue = 240;
    const saturation = 70;
    const lightness = 100 - (scaled * 60);
    return {
        bg: `hsl(${hue}, ${saturation}%, ${lightness}%)`,
        text: lightness < 60 ? 'white' : 'black'
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
    heatmapData.taxa.forEach(t => {
        const taxonName = t[rank] || 'Unknown';
        if (!aggregated[taxonName]) {
            aggregated[taxonName] = new Array(heatmapData.samples.length).fill(0);
        }
        t.abundances.forEach((val, i) => {
            aggregated[taxonName][i] += val;
        });
    });

    const allTaxa = Object.keys(aggregated).map(name => ({
        name: name,
        sum: sampleIndices.reduce((acc, idx) => acc + (aggregated[name][idx] || 0), 0),
        abundances: aggregated[name]
    }))
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
        html += `<tr><td class="p-2 font-medium ${isOther ? 'bg-gray-50 italic' : 'bg-white'} border sticky left-0 z-10 whitespace-nowrap shadow-sm">${t.name}</td>`;
        sampleIndices.forEach(idx => {
            const val = t.abundances[idx] || 0;
            const { bg, text } = getHeatmapColor(val);
            html += `<td class="p-2 text-center border heatmap-cell" style="background-color: ${bg}; color: ${text}" title="${t.name} in ${heatmapData.samples[idx]}: ${(val * 100).toFixed(4)}%">
                        ${val > 0.001 ? (val * 100).toFixed(1) + '%' : (val > 0 ? '<0.1%' : '-')}
                     </td>`;
        });
        html += '</tr>';
    });
    html += '</tbody></table>';
    container.innerHTML = html;
}

// --- Chart Logic ---
function getColor(name, index) {
    if (name.startsWith('Other')) return '#9ca3af';
    const colors = ['#4f46e5', '#ef4444', '#10b981', '#f59e0b', '#3b82f6', '#8b5cf6', '#ec4899', '#6366f1', '#14b8a6', '#f97316'];
    return colors[index % colors.length];
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
            const taxon = d[rank] || 'Unknown';
            globalTotals[taxon] = (globalTotals[taxon] || 0) + d.abundance;
            if (!presence[taxon]) presence[taxon] = new Set();
            if (d.abundance > 0) presence[taxon].add(s.name);
        });
    });

    const commonTaxa = Object.keys(globalTotals).filter(t => presence[t].size === selectedSamples.length);
    const topTaxa = commonTaxa.sort((a, b) => globalTotals[b] - globalTotals[a]).slice(0, topN);

    const datasets = topTaxa.map((taxon, i) => {
        const taxonData = selectedSamples.map(s => {
            const match = s.data.filter(d => d[rank] === taxon);
            return { abundance: match.reduce((sum, d) => sum + d.abundance, 0) * 100, counts: Math.round(match.reduce((sum, d) => sum + d.counts, 0)) };
        });
        return { label: taxon, data: taxonData.map(d => d.abundance), counts: taxonData.map(d => d.counts), backgroundColor: getColor(taxon, i) };
    });

    const otherData = selectedSamples.map(s => {
        const other = s.data.filter(d => !topTaxa.includes(d[rank]));
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

// Initialize
populateSampleSelector();
updateChart();
renderHeatmap();

// Initial sort for summary table
setTimeout(() => sortSummaryTable(0), 100); 
