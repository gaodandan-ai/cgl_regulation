(function (global) {
    const OBJECTIVE_CATALOG = {
        biomass: {
            label: "Growth / biomass objective",
            reactionId: null,
            type: "model_default",
            description: "Use the model's default biomass objective."
        },
        glutamate_export: {
            label: "Glutamate export objective",
            reactionId: null,
            type: "requires_reaction_selection",
            description: "Requires a verified glutamate exchange/export reaction."
        },
        lysine_export: {
            label: "Lysine export objective",
            reactionId: null,
            type: "requires_reaction_selection",
            description: "Requires a verified lysine exchange/export reaction."
        }
    };

    global.modelObjectives = {
        OBJECTIVE_CATALOG
    };
})(window);
