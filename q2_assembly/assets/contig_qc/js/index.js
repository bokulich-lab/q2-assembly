$(document).ready(function () {
    removeBS3refs();
    adjustTagsToBS3();

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
    const initialCategory = categories[0];
    updateValueDropdown(initialCategory, 'All');
    document.getElementById('categoryDropdown').value = initialCategory;

    // Render each plot in its own container
    let vegaViews = {};
    vegaEmbed('#vega-contig-length', vegaContigLengthSpec, {actions: true}).then(res => {
        vegaViews.contigLength = res.view;
        vegaViews.contigLength.signal('category_param', initialCategory)
            .signal('value_param', 'All')
            .runAsync();
    });
    vegaEmbed('#vega-nx-curve', vegaNxCurveSpec, {actions: true}).then(res => {
        vegaViews.nxCurve = res.view;
        vegaViews.nxCurve.signal('category_param', initialCategory)
            .signal('value_param', 'All')
            .runAsync();
    });
    vegaEmbed('#vega-gc-content', vegaGcContentSpec, {actions: true}).then(res => {
        vegaViews.gcContent = res.view;
        vegaViews.gcContent.signal('category_param', initialCategory)
            .signal('value_param', 'All')
            .runAsync();
    });
    vegaEmbed('#vega-cumulative-length', vegaCumulativeLengthSpec, {actions: true}).then(res => {
        vegaViews.cumulativeLength = res.view;
        vegaViews.cumulativeLength.signal('category_param', initialCategory)
            .signal('value_param', 'All')
            .runAsync();
    });

    // Remove the loading spinner
    document.getElementById('loading').remove();

    // Dropdown event listeners
    document.getElementById('categoryDropdown').addEventListener('change', function () {
        const category = this.value;
        updateValueDropdown(category, 'All');
        // Update Vega params for all plots
        Object.values(vegaViews).forEach(view => {
            if (view) {
                view.signal('category_param', category).signal('value_param', 'All').runAsync();
            }
        });
    });
    document.getElementById('valueDropdown').addEventListener('change', function () {
        const value = this.value;
        Object.values(vegaViews).forEach(view => {
            if (view) {
                view.signal('value_param', value).runAsync();
            }
        });
    });

});