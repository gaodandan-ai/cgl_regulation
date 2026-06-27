import os
import logging
import cobra

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("model_loader")

MODEL_PATH = os.path.join(os.path.dirname(__file__), "models", "iCW773.xml")

# Cached model instance and status
_cached_model = None
_load_error = None
_initialized = False

def get_model_status():
    """
    Get the status of the model file and if it has loaded successfully.
    """
    global _cached_model, _load_error, _initialized
    
    # Check if file exists
    if not os.path.exists(MODEL_PATH):
        return {
            "loaded": False,
            "error": f"Model file missing at backend/models/iCW773.xml",
            "model_id": None,
            "reaction_count": 0,
            "gene_count": 0,
            "metabolite_count": 0
        }
        
    if _load_error:
        return {
            "loaded": False,
            "error": _load_error,
            "model_id": None,
            "reaction_count": 0,
            "gene_count": 0,
            "metabolite_count": 0
        }
        
    if _cached_model:
        return {
            "loaded": True,
            "error": None,
            "model_id": _cached_model.id,
            "reaction_count": len(_cached_model.reactions),
            "gene_count": len(_cached_model.genes),
            "metabolite_count": len(_cached_model.metabolites)
        }
        
    # Attempt lazy loading on status check if not loaded yet
    try:
        load_model_if_needed()
        return {
            "loaded": True,
            "error": None,
            "model_id": _cached_model.id,
            "reaction_count": len(_cached_model.reactions),
            "gene_count": len(_cached_model.genes),
            "metabolite_count": len(_cached_model.metabolites)
        }
    except Exception as e:
        return {
            "loaded": False,
            "error": str(e),
            "model_id": None,
            "reaction_count": 0,
            "gene_count": 0,
            "metabolite_count": 0
        }

def load_model_if_needed():
    """
    Load the model from disk if it hasn't been loaded already.
    Returns the model or raises an Exception if loading fails.
    """
    global _cached_model, _load_error, _initialized
    
    if _initialized:
        if _load_error:
            raise Exception(_load_error)
        return _cached_model

    _initialized = True
    
    if not os.path.exists(MODEL_PATH):
        _load_error = f"Model file missing at backend/models/iCW773.xml"
        logger.error(_load_error)
        raise Exception(_load_error)
        
    try:
        logger.info(f"Loading SBML model from {MODEL_PATH}...")
        _cached_model = cobra.io.read_sbml_model(MODEL_PATH)
        logger.info("Model loaded successfully!")
        _load_error = None
        return _cached_model
    except Exception as e:
        _load_error = f"Failed to parse SBML model: {str(e)}"
        logger.error(_load_error)
        raise Exception(_load_error)
