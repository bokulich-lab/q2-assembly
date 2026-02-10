$(document).ready(function () {
    removeBS3refs();
    adjustTagsToBS3();

    hasClientMetadata = Boolean(hasClientMetadata)

    // hasClientMetadata, sampleIdsByMetadata, categories, values, vega specs are injected from Python by the Jinja template
    let vegaViews = {};

    const customLegendDiv = document.getElementById('custom-legend');
    const selectAllBtn = document.getElementById('selectAllLegendBtn');
    const unselectAllBtn = document.getElementById('unselectAllLegendBtn');

    // These will be assigned if hasClientMetadata is true and elements are found
    let categoryDropdown = null;
    let valueDropdown = null;

    function debounce(func, delay) {
        let timeout;
        return function(...args) {
            const context = this;
            clearTimeout(timeout);
            timeout = setTimeout(() => func.apply(context, args), delay);
        };
    }

    function refreshVegaViews() {
        let currentCategorySignal = "None"; // Default for Vega signal if no metadata/dropdowns
        let currentValueSignal = "All";    // Default for Vega signal if no metadata/dropdowns

        if (hasClientMetadata && categoryDropdown && valueDropdown) { // Check if dropdowns were initialized
            currentCategorySignal = categoryDropdown.value;
            currentValueSignal = valueDropdown.value;
        }

        const selectedLegendItems = document.querySelectorAll('#custom-legend .custom-legend-item.selected');
        const selectedSamples = Array.from(selectedLegendItems).map(item => {
            return { sample: item.dataset.value };
        });

        Object.values(vegaViews).forEach(view => {
            if (view) {
                view.signal('category_param', currentCategorySignal)
                    .signal('value_param', currentValueSignal)
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

        customLegendDiv.textContent = '';
        const fragment = document.createDocumentFragment();
        let legendKeys = [];

        if (selectedValue === 'All') {
            if (selectedCategory && sampleIdsByMetadata[selectedCategory]) {
                let samplesForCategorySet = new Set();
                for (const valueKey in sampleIdsByMetadata[selectedCategory]) {
                    if (Object.prototype.hasOwnProperty.call(sampleIdsByMetadata[selectedCategory], valueKey)) {
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
                legendKeys = sampleIdsByMetadata.all_samples;
            }
        } else {
            if (selectedCategory && sampleIdsByMetadata[selectedCategory] &&
                Object.prototype.hasOwnProperty.call(sampleIdsByMetadata[selectedCategory], selectedValue)) {
                legendKeys = sampleIdsByMetadata[selectedCategory][selectedValue] || [];
            }
        }

        if (legendKeys.length === 0) {
            return;
        }

        const n = legendKeys.length;
        const colorRange = n > 0
          ? d3.quantize(d3.interpolateViridis, n + 2).slice(1, -1)
          : [];

        const viridisColorScale = d3.scaleOrdinal()
          .domain(legendKeys)
          .range(colorRange);

        legendKeys.forEach((key) => {
            const legendItem = document.createElement('span');
            legendItem.classList.add('legend-item', 'custom-legend-item', 'selected');
            legendItem.dataset.value = key;

            const colorSwatch = document.createElement('span');
            colorSwatch.classList.add('legend-color-swatch');
            colorSwatch.style.backgroundColor = viridisColorScale(key);

            const textSpan = document.createElement('span');
            textSpan.textContent = key;

            legendItem.appendChild(colorSwatch);
            legendItem.appendChild(textSpan);
            fragment.appendChild(legendItem);
        });
        customLegendDiv.appendChild(fragment);
    }

    // This helper populates the value dropdown; only call if valueDropdown exists.
    function populateValueDropdownWithOptions(category, selectedValue) {
        if (!valueDropdown) return; // Should only be called if valueDropdown element is available

        valueDropdown.textContent = '';
        const allOpt = document.createElement('option');
        allOpt.value = 'All';
        allOpt.textContent = 'All';
        if (selectedValue === 'All') allOpt.selected = true;
        valueDropdown.appendChild(allOpt);

        // 'values' is injected from Python
        if (category && values[category]) {
            values[category].forEach(v => {
                const opt = document.createElement('option');
                opt.value = v;
                opt.textContent = v;
                if (v === selectedValue) opt.selected = true;
                valueDropdown.appendChild(opt);
            });
        }
    }

    if (hasClientMetadata) {
        categoryDropdown = document.getElementById('categoryDropdown');
        valueDropdown = document.getElementById('valueDropdown');

        if (categoryDropdown && valueDropdown) { // Ensure dropdowns were actually found in DOM
            // Populate category dropdown ('categories' is injected from Python)
            categories.forEach(cat => {
                const opt = document.createElement('option');
                opt.value = cat;
                opt.textContent = cat.charAt(0).toUpperCase() + cat.slice(1);
                categoryDropdown.appendChild(opt);
            });

            const initialCategory = categories.length > 0 ? categories[0] : null;
            if (initialCategory) {
                categoryDropdown.value = initialCategory;
                populateValueDropdownWithOptions(initialCategory, 'All');
            } else {
                // hasClientMetadata is true, but no categories (e.g., metadata object exists but is empty)
                populateValueDropdownWithOptions(null, 'All'); // Value dropdown will show only 'All'
            }

            // Event listeners for dropdowns
            categoryDropdown.addEventListener('change', function () {
                const category = this.value;
                populateValueDropdownWithOptions(category, 'All');
                updateCustomLegend(category, 'All');
                debouncedRefreshVegaViews();
            });

            valueDropdown.addEventListener('change', function () {
                const value = this.value;
                const currentCategory = categoryDropdown.value; // categoryDropdown is checked to exist here
                updateCustomLegend(currentCategory, value);
                debouncedRefreshVegaViews();
            });
        } else {
            // This case means index.html indicated metadata, but dropdowns weren't found.
            // Could happen if HTML structure is changed without updating JS, or IDs are wrong.
            console.warn("Metadata dropdowns (categoryDropdown/valueDropdown) not found in DOM even though hasClientMetadata is true.");
        }
    }

    // Initial legend update:
    // If metadata dropdowns were set up, use the current category from categoryDropdown.
    // Otherwise (no metadata or dropdowns not found), selectedCategory will be null.
    const initialLegendCategory = (hasClientMetadata && categoryDropdown) ? categoryDropdown.value : null;
    updateCustomLegend(initialLegendCategory, 'All');

    const vegaEmbedPromises = [];
    vegaEmbedPromises.push(
        vegaEmbed('#vega-contig-length', vegaContigLengthSpec, {actions: true, renderer: 'canvas'}).then(res => {
            vegaViews.contigLength = res.view;
            document.getElementById('spinner-contig-length')?.remove();
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

    Promise.all(vegaEmbedPromises).then(() => {
        refreshVegaViews();
    }).catch(error => {
        console.error("Error embedding Vega views:", error);
        refreshVegaViews();
    });

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

    if (customLegendDiv) {
        customLegendDiv.addEventListener('click', function(event) {
            const targetItem = event.target.closest('.custom-legend-item');
            if (!targetItem) return;
            targetItem.classList.toggle('selected');
            debouncedRefreshVegaViews();
        });
    }
});