$(document).ready(function () {
    removeBS3refs();
    adjustTagsToBS3();

    // sampleIdsByMetadata, categories, values, vega specs are injected from Python by the Jinja template
    let vegaViews = {}; // Initialize vegaViews here

    function updateVegaHighlightSignal() {
        const selectedLegendItems = document.querySelectorAll('#custom-legend .custom-legend-item.selected');
        const selectedSamples = Array.from(selectedLegendItems).map(item => {
            return { sample: item.dataset.value }; // Format for Vega multi-selection
        });

        Object.values(vegaViews).forEach(view => {
            if (view) {
                view.signal('highlight_samples_param', selectedSamples).runAsync();
            }
        });
    }

    function updateCustomLegend(selectedCategory, selectedValue) {
        const customLegendDiv = document.getElementById('custom-legend');
        if (!customLegendDiv || typeof sampleIdsByMetadata === 'undefined') {
            if (customLegendDiv) customLegendDiv.innerHTML = 'Error: sampleIdsByMetadata not loaded.';
            return;
        }

        customLegendDiv.innerHTML = ''; // Clear previous legend items
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
            // If no specific category led to keys, or if category itself was not found/null,
            // or if legendKeys is still empty after processing a category (e.g. category has no samples)
            // default to all samples if available.
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

        // Define the color scale based on the final legendKeys
        // Use d3.quantize with d3.interpolateViridis to get a dynamic range
        // of colors spread across the Viridis spectrum, based on the number of legend keys.
        // const colorRange = legendKeys.length > 0 ? d3.quantize(d3.interpolateViridis, legendKeys.length) : [];
        // const viridisColorScale = d3.scaleOrdinal().domain(legendKeys).range(colorRange);

        const n = legendKeys.length;
        // quantize at n+2 points, then drop the first and last
        const colorRange = n > 0
          ? d3.quantize(d3.interpolateViridis, n + 2).slice(1, -1)
          : [];

        const viridisColorScale = d3.scaleOrdinal()
          .domain(legendKeys)
          .range(colorRange);

        console.log(viridisColorScale(legendKeys[0]));
        console.log(viridisColorScale(legendKeys[legendKeys.length - 1]))
        legendKeys.forEach((key) => { // key is now a sample ID
            const legendItem = document.createElement('span');
            legendItem.classList.add('legend-item', 'custom-legend-item', 'selected'); // Initially selected
            legendItem.style.marginRight = '15px';
            legendItem.style.marginBottom = '5px';
            legendItem.style.display = 'flex';
            legendItem.style.alignItems = 'center';
            legendItem.dataset.value = key; // Store the sample ID

            const colorSwatch = document.createElement('span');
            colorSwatch.style.display = 'inline-block';
            colorSwatch.style.width = '12px';
            colorSwatch.style.height = '12px';
            // Use the D3 scale function to get the color
            colorSwatch.style.backgroundColor = viridisColorScale(key);
            colorSwatch.style.borderRadius = '2px';
            colorSwatch.style.marginRight = '5px';

            const textSpan = document.createElement('span');
            textSpan.textContent = key; // Display sample ID

            legendItem.appendChild(colorSwatch);
            legendItem.appendChild(textSpan);
            customLegendDiv.appendChild(legendItem);
        });
    }

    // Populate category dropdown
    const categoryDropdown = document.getElementById('categoryDropdown');
    categories.forEach(cat => {
        const opt = document.createElement('option');
        opt.value = cat;
        opt.textContent = cat.charAt(0).toUpperCase() + cat.slice(1);
        categoryDropdown.appendChild(opt);
    });

    // Helper to update value dropdown
    function updateValueDropdown(category, selectedValue) {
        const valueDropdown = document.getElementById('valueDropdown');
        valueDropdown.innerHTML = '';
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
        // updateVegaHighlightSignal(); // Defer this call
    } else {
        // No categories from metadata. Value dropdown will be empty except for "All".
        updateValueDropdown(null, 'All');
        updateCustomLegend(null, 'All'); // Should fall back to sampleIdsByMetadata.all_samples
        // updateVegaHighlightSignal(); // Defer this call
    }

    // Render each plot in its own container
    vegaEmbedPromises.push(
        vegaEmbed('#vega-contig-length', vegaContigLengthSpec, {actions: true}).then(res => {
            vegaViews.contigLength = res.view;
            document.getElementById('spinner-contig-length')?.remove();
            if (initialCategory) {
                vegaViews.contigLength.signal('category_param', initialCategory)
                    .signal('value_param', 'All')
                    .runAsync();
            } else {
                 vegaViews.contigLength.signal('value_param', 'All').runAsync();
            }
        })
    );
    vegaEmbedPromises.push(
        vegaEmbed('#vega-nx-curve', vegaNxCurveSpec, {actions: true}).then(res => {
            vegaViews.nxCurve = res.view;
            document.getElementById('spinner-nx-curve')?.remove();
            if (initialCategory) {
                vegaViews.nxCurve.signal('category_param', initialCategory)
                    .signal('value_param', 'All')
                    .runAsync();
            } else {
                vegaViews.nxCurve.signal('value_param', 'All').runAsync();
            }
        })
    );
    vegaEmbedPromises.push(
        vegaEmbed('#vega-gc-content', vegaGcContentSpec, {actions: true}).then(res => {
            vegaViews.gcContent = res.view;
            document.getElementById('spinner-gc-content')?.remove();
            if (initialCategory) {
                vegaViews.gcContent.signal('category_param', initialCategory)
                    .signal('value_param', 'All')
                    .runAsync();
            } else {
                vegaViews.gcContent.signal('value_param', 'All').runAsync();
            }
        })
    );
    vegaEmbedPromises.push(
        vegaEmbed('#vega-cumulative-length', vegaCumulativeLengthSpec, {actions: true}).then(res => {
            vegaViews.cumulativeLength = res.view;
            document.getElementById('spinner-cumulative-length')?.remove();
            if (initialCategory) {
                vegaViews.cumulativeLength.signal('category_param', initialCategory)
                    .signal('value_param', 'All')
                    .runAsync();
            } else {
                vegaViews.cumulativeLength.signal('value_param', 'All').runAsync();
            }
        })
    );

    // Global spinner removal is no longer needed here
    // document.getElementById('loading').remove(); 

    // Wait for all Vega plots to be embedded before the first highlight signal update
    Promise.all(vegaEmbedPromises).then(() => {
        updateVegaHighlightSignal(); // Set initial highlight state now
    }).catch(error => {
        console.error("Error embedding Vega views:", error);
        // Optionally, still call updateVegaHighlightSignal or handle error state
        updateVegaHighlightSignal(); // Attempt to update with whatever views are available
    });

    // Event listeners for Select/Unselect All buttons for the custom legend
    const selectAllBtn = document.getElementById('selectAllLegendBtn');
    const unselectAllBtn = document.getElementById('unselectAllLegendBtn');

    if (selectAllBtn) {
        selectAllBtn.addEventListener('click', function() {
            const legendItems = document.querySelectorAll('#custom-legend .custom-legend-item');
            legendItems.forEach(item => {
                item.classList.add('selected');
            });
            updateVegaHighlightSignal();
        });
    }

    if (unselectAllBtn) {
        unselectAllBtn.addEventListener('click', function() {
            const legendItems = document.querySelectorAll('#custom-legend .custom-legend-item');
            legendItems.forEach(item => {
                item.classList.remove('selected');
            });
            updateVegaHighlightSignal();
        });
    }

    // Dropdown event listeners
    document.getElementById('categoryDropdown').addEventListener('change', function () {
        const category = this.value;
        updateValueDropdown(category, 'All');
        updateCustomLegend(category, 'All'); 
        updateVegaHighlightSignal(); // Ensure plots reflect new legend's all-selected state
        // Update Vega params for all plots
        Object.values(vegaViews).forEach(view => {
            if (view) {
                view.signal('category_param', category).signal('value_param', 'All').runAsync();
            }
        });
    });
    document.getElementById('valueDropdown').addEventListener('change', function () {
        const value = this.value;
        const currentCategory = document.getElementById('categoryDropdown').value;
        updateCustomLegend(currentCategory, value);
        updateVegaHighlightSignal(); // Ensure plots reflect new legend's all-selected state
        Object.values(vegaViews).forEach(view => {
            if (view) {
                view.signal('value_param', value).runAsync();
            }
        });
    });

    // Custom Legend Interactivity
    const customLegend = document.getElementById('custom-legend');
    if (customLegend) {
        const legendItems = customLegend.querySelectorAll('.custom-legend-item');

        // Initially, mark all legend items as selected
        legendItems.forEach(item => {
            item.classList.add('selected');
        });

        customLegend.addEventListener('click', function(event) {
            const targetItem = event.target.closest('.custom-legend-item');
            if (targetItem) {
                targetItem.classList.toggle('selected');
                updateVegaHighlightSignal(); // Update Vega plots based on new selection
            }
        });
    }

});