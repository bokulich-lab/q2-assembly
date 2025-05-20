$(document).ready(function () {
    removeBS3refs();
    adjustTagsToBS3();

    // sampleIdsByMetadata, categories, values, vega specs are injected from Python by the Jinja template
    let vegaViews = {}; // Initialize vegaViews here

    // Cache frequently used DOM nodes
    const customLegendDiv = document.getElementById('custom-legend');
    const categoryDropdown = document.getElementById('categoryDropdown');
    const valueDropdown = document.getElementById('valueDropdown');
    const selectAllBtn = document.getElementById('selectAllLegendBtn');
    const unselectAllBtn = document.getElementById('unselectAllLegendBtn');

    // Debounce function
    function debounce(func, delay) {
        let timeout;
        return function(...args) {
            const context = this;
            clearTimeout(timeout);
            timeout = setTimeout(() => func.apply(context, args), delay);
        };
    }

    // Central function to update all Vega views via signals
    function refreshVegaViews() {
        const currentCategory = categoryDropdown.value;
        const currentValue = valueDropdown.value;

        const selectedLegendItems = document.querySelectorAll('#custom-legend .custom-legend-item.selected');
        const selectedSamples = Array.from(selectedLegendItems).map(item => {
            return { sample: item.dataset.value }; // Format for Vega multi-selection
        });

        // Optional: Log current signal values for debugging
        // console.log('Refreshing Vega Views. Signals:', {
        //     category: currentCategory,
        //     value: currentValue,
        //     highlighted: selectedSamples.map(s => s.sample)
        // });

        Object.values(vegaViews).forEach(view => {
            if (view) {
                view.signal('category_param', currentCategory)
                    .signal('value_param', currentValue)
                    .signal('highlight_samples_param', selectedSamples)
                    .runAsync();
            }
        });
    }

    const debouncedRefreshVegaViews = debounce(refreshVegaViews, 150);

    function updateCustomLegend(selectedCategory, selectedValue) {
        if (!customLegendDiv || typeof sampleIdsByMetadata === 'undefined') {
            if (customLegendDiv) customLegendDiv.innerHTML = 'Error: sampleIdsByMetadata not loaded.';
            return;
        }

        customLegendDiv.textContent = ''; // Clear previous legend items
        const fragment = document.createDocumentFragment();
        let legendKeys = [];

        if (selectedValue === 'All') {
            if (selectedCategory && sampleIdsByMetadata[selectedCategory]) {
                // Collect all unique sample IDs under the specific category's values
                let samplesForCategorySet = new Set();
                for (const valueKey in sampleIdsByMetadata[selectedCategory]) {
                    if (Object.prototype.hasOwnProperty.call(sampleIdsByMetadata[selectedCategory], valueKey)) {
                        // Ensure sampleIdsByMetadata[selectedCategory][valueKey] is an array
                        if (Array.isArray(sampleIdsByMetadata[selectedCategory][valueKey])) {
                            sampleIdsByMetadata[selectedCategory][valueKey].forEach(sampleId => {
                                samplesForCategorySet.add(sampleId);
                            });
                        }
                    }
                }
                legendKeys = Array.from(samplesForCategorySet).sort();
            }
            if (legendKeys.length === 0 && sampleIdsByMetadata.all_samples && sampleIdsByMetadata.all_samples.length > 0) {
                legendKeys = sampleIdsByMetadata.all_samples; // Already sorted from Python
            }
        } else {
            // Specific category and value selected
            if (selectedCategory && sampleIdsByMetadata[selectedCategory] &&
                Object.prototype.hasOwnProperty.call(sampleIdsByMetadata[selectedCategory], selectedValue)) {
                legendKeys = sampleIdsByMetadata[selectedCategory][selectedValue] || []; // Already sorted from Python
            }
        }

        if (legendKeys.length === 0) {
            // customLegendDiv.textContent = 'No legend items to display.';
            return;
        }

        const n = legendKeys.length;
        // quantize at n+2 points, then drop the first and last
        const colorRange = n > 0
          ? d3.quantize(d3.interpolateViridis, n + 2).slice(1, -1)
          : [];

        const viridisColorScale = d3.scaleOrdinal()
          .domain(legendKeys)
          .range(colorRange);

        legendKeys.forEach((key) => { // key is now a sample ID
            const legendItem = document.createElement('span');
            legendItem.classList.add('legend-item', 'custom-legend-item', 'selected'); // Initially selected
            legendItem.dataset.value = key; // Store the sample ID

            const colorSwatch = document.createElement('span');
            colorSwatch.classList.add('legend-color-swatch');
            colorSwatch.style.backgroundColor = viridisColorScale(key);

            const textSpan = document.createElement('span');
            textSpan.textContent = key; // Display sample ID

            legendItem.appendChild(colorSwatch);
            legendItem.appendChild(textSpan);
            fragment.appendChild(legendItem);
        });
        customLegendDiv.appendChild(fragment);
    }

    // Populate category dropdown
    categories.forEach(cat => {
        const opt = document.createElement('option');
        opt.value = cat;
        opt.textContent = cat.charAt(0).toUpperCase() + cat.slice(1);
        categoryDropdown.appendChild(opt);
    });

    // Helper to update value dropdown
    function updateValueDropdown(category, selectedValue) {
        valueDropdown.textContent = '';
        // Add 'All' option
        const allOpt = document.createElement('option');
        allOpt.value = 'All';
        allOpt.textContent = 'All';
        if (selectedValue === 'All') allOpt.selected = true;
        valueDropdown.appendChild(allOpt);
        if (values[category]) {
            values[category].forEach(v => {
                const opt = document.createElement('option');
                opt.value = v;
                opt.textContent = v;
                if (v === selectedValue) opt.selected = true;
                valueDropdown.appendChild(opt);
            });
        }
    }

    // Initial dropdown population
    const initialCategory = categories.length > 0 ? categories[0] : null;

    // Store all vegaEmbed promises
    const vegaEmbedPromises = [];

    if (initialCategory) {
        updateValueDropdown(initialCategory, 'All');
        document.getElementById('categoryDropdown').value = initialCategory;
        updateCustomLegend(initialCategory, 'All'); // Generate initial legend
        // refreshVegaViews(); // Called after all plots are embedded
    } else {
        // No categories from metadata. Value dropdown will be empty except for "All".
        updateValueDropdown(null, 'All');
        updateCustomLegend(null, 'All'); // Should fall back to sampleIdsByMetadata.all_samples
        // refreshVegaViews(); // Called after all plots are embedded
    }

    // Render each plot in its own container
    vegaEmbedPromises.push(
        vegaEmbed('#vega-contig-length', vegaContigLengthSpec, {actions: true, renderer: 'canvas'}).then(res => {
            vegaViews.contigLength = res.view;
            document.getElementById('spinner-contig-length')?.remove();
            // Initial signal setting will be handled by refreshVegaViews after all embeds
        })
    );
    vegaEmbedPromises.push(
        vegaEmbed('#vega-nx-curve', vegaNxCurveSpec, {actions: true, renderer: 'canvas'}).then(res => {
            vegaViews.nxCurve = res.view;
            document.getElementById('spinner-nx-curve')?.remove();
        })
    );
    vegaEmbedPromises.push(
        vegaEmbed('#vega-gc-content', vegaGcContentSpec, {actions: true, renderer: 'canvas'}).then(res => {
            vegaViews.gcContent = res.view;
            document.getElementById('spinner-gc-content')?.remove();
        })
    );
    vegaEmbedPromises.push(
        vegaEmbed('#vega-cumulative-length', vegaCumulativeLengthSpec, {actions: true, renderer: 'canvas'}).then(res => {
            vegaViews.cumulativeLength = res.view;
            document.getElementById('spinner-cumulative-length')?.remove();
        })
    );

    // Global spinner removal is no longer needed here
    // document.getElementById('loading').remove();

    // Wait for all Vega plots to be embedded before the first signal update
    Promise.all(vegaEmbedPromises).then(() => {
        refreshVegaViews(); // Set initial signal state for all charts
    }).catch(error => {
        console.error("Error embedding Vega views:", error);
        // Optionally, still attempt to set signals or handle error state
        refreshVegaViews(); // Attempt to update with whatever views are available
    });

    // Event listeners for Select/Unselect All buttons for the custom legend
    if (selectAllBtn) {
        selectAllBtn.addEventListener('click', function() {
            const legendItems = customLegendDiv.querySelectorAll('.custom-legend-item');
            legendItems.forEach(item => {
                item.classList.add('selected');
            });
            debouncedRefreshVegaViews();
        });
    }

    if (unselectAllBtn) {
        unselectAllBtn.addEventListener('click', function() {
            const legendItems = customLegendDiv.querySelectorAll('.custom-legend-item');
            legendItems.forEach(item => {
                item.classList.remove('selected');
            });
            debouncedRefreshVegaViews();
        });
    }

    // Dropdown event listeners
    document.getElementById('categoryDropdown').addEventListener('change', function () {
        const category = this.value;
        updateValueDropdown(category, 'All'); // This will set valueDropdown to 'All'
        updateCustomLegend(category, 'All'); // Rebuilds legend, all new items are selected by default
        debouncedRefreshVegaViews();
    });
    document.getElementById('valueDropdown').addEventListener('change', function () {
        const value = this.value;
        const currentCategory = document.getElementById('categoryDropdown').value;
        updateCustomLegend(currentCategory, value); // Rebuilds legend, all new items are selected by default
        debouncedRefreshVegaViews();
    });

    // Custom Legend Interactivity (event delegation)
    if (customLegendDiv) {
        customLegendDiv.addEventListener('click', function(event) {
            const targetItem = event.target.closest('.custom-legend-item');
            if (!targetItem) return;
            targetItem.classList.toggle('selected');
            debouncedRefreshVegaViews();
        });
    }

});