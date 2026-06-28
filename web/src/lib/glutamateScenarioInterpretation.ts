export function generateGlutamateScenarioInterpretation(result: any, selectedReaction: any): string {
  const classification = selectedReaction?.classification || "uncertain";
  const reactionId = selectedReaction?.reactionId || "";
  
  let labelPhrase = "selected glutamate-associated reaction";
  if (classification === "exchange" || classification === "export") {
    labelPhrase = "selected glutamate export/exchange reaction";
  }
  
  let fluxChange = 0;
  if (result.trackedFluxes && Array.isArray(result.trackedFluxes)) {
    const tf = result.trackedFluxes.find((f: any) => f.reactionId === reactionId);
    if (tf) {
      fluxChange = tf.fluxChange || 0;
    }
  }
  
  let effectPhrase = "";
  if (fluxChange > 1e-5) {
    effectPhrase = `The ${labelPhrase} (${reactionId}) shows increased predicted flux under this perturbation.`;
  } else if (fluxChange < -1e-5) {
    effectPhrase = `The ${labelPhrase} (${reactionId}) shows reduced predicted flux under this perturbation.`;
  } else {
    effectPhrase = `The ${labelPhrase} (${reactionId}) shows no significant predicted flux change under this perturbation.`;
  }
  
  return `${effectPhrase} This is an in silico prediction and requires experimental validation.`;
}
