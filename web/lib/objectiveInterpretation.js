(function (global) {
    function generateObjectiveSimulationInterpretation(result) {
        const objectiveType = result.objective?.objectiveType || 'biomass';
        
        let intro = '';
        if (objectiveType === 'biomass') {
            intro = 'under the current biomass objective';
        } else {
            intro = 'under the selected reaction objective';
        }
        
        let glutamateIncreased = false;
        let trackedGluReaction = '';
        
        if (result.trackedFluxes && Array.isArray(result.trackedFluxes)) {
            for (const tf of result.trackedFluxes) {
                const idLower = (tf.reactionId || '').toLowerCase();
                if (idLower.includes('glu')) {
                    trackedGluReaction = tf.reactionId;
                    if (tf.perturbedFlux > tf.baselineFlux + 0.0001) {
                        glutamateIncreased = true;
                    }
                }
            }
        }
        
        let body = 'this perturbation is predicted to alter downstream metabolic flux.';
        if (glutamateIncreased && trackedGluReaction) {
            body = `the selected glutamate-associated reaction (${trackedGluReaction}) shows increased predicted flux.`;
        } else if (trackedGluReaction) {
            body = 'this perturbation is predicted to alter glutamate-associated flux.';
        }
        
        return `Under the selected model constraints and ${intro}, ${body} This is an in silico prediction and requires experimental validation.`;
    }

    global.objectiveInterpretation = {
        generateObjectiveSimulationInterpretation
    };
})(window);
