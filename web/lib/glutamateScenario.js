(function (global) {
    function classifyGlutamateReaction(reaction) {
        const rxnId = (reaction.reactionId || reaction.id || '').toLowerCase();
        const name = (reaction.name || '').toLowerCase();
        const equation = (reaction.equation || reaction.reaction || '').toLowerCase();
        
        const hasExtracellular = equation.includes('glu__l_e') || equation.includes('glu_e');
        const hasIntracellular = equation.includes('glu__l_c') || equation.includes('glu_c');
        
        const isExchange = rxnId.startsWith('ex_') || rxnId.includes('_ex');
        
        if (isExchange) {
            return {
                classification: 'exchange',
                confidence: 'high',
                reason: 'Reaction ID suggests exchange and equation represents extracellular L-glutamate boundary flux.'
            };
        } else if (name.includes('export') || rxnId.includes('export') || name.includes('secretion')) {
            return {
                classification: 'export',
                confidence: 'high',
                reason: 'Reaction name or equation explicitly suggests extracellular glutamate secretion or export.'
            };
        } else if (name.includes('transport') || (hasExtracellular && hasIntracellular)) {
            return {
                classification: 'transport',
                confidence: 'medium',
                reason: 'Reaction represents transport of L-glutamate across cellular compartments.'
            };
        } else if (hasIntracellular && !hasExtracellular) {
            if (name.includes('synth') || name.includes('dehydrogenase') || name.includes('transaminase')) {
                return {
                    classification: 'biosynthesis',
                    confidence: 'medium',
                    reason: 'Intracellular enzymatic reaction converting reactants to L-glutamate.'
                };
            } else if (name.includes('decarboxylase') || name.includes('kinase') || name.includes('synthase')) {
                return {
                    classification: 'consumption',
                    confidence: 'medium',
                    reason: 'Intracellular reaction consuming L-glutamate.'
                };
            } else {
                return {
                    classification: 'uncertain',
                    confidence: 'low',
                    reason: 'Intracellular glutamate conversion reaction of uncertain direction.'
                };
            }
        } else {
            return {
                classification: 'uncertain',
                confidence: 'low',
                reason: 'Glutamate-associated reaction of uncertain category or compartment.'
            };
        }
    }

    // Default global state
    const glutamateState = {
        selectedGlutamateReactionId: null,
        selectedGlutamateReactionClass: null,
        userVerified: false
    };

    global.glutamateScenario = {
        classifyGlutamateReaction,
        glutamateState
    };
})(window);
