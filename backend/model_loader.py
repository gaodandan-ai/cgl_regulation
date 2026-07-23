import os
import logging
import threading

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("model_loader")

# Thermodynamic pruning (applied once after model load)
try:
    from thermo_pruner import apply_thermodynamic_pruning, get_pruning_report
    _THERMO_AVAILABLE = True
except ImportError:
    _THERMO_AVAILABLE = False
    logger.warning("thermo_pruner not found – model loaded without thermodynamic pruning.")

MODEL_PATH = os.path.join(os.path.dirname(__file__), "models", "iCW773.xml")

# Cached model instance and status
_cached_model = None
_load_error = None
_initialized = False
_loading = False
_model_condition = threading.Condition()
_preload_thread = None

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
    global _cached_model, _load_error, _initialized, _loading

    # Model loading can be started by the background warm-up thread while a
    # simulation request arrives. Make the first load single-flight and let
    # callers wait for that same result instead of parsing the SBML twice.
    with _model_condition:
        while _loading:
            _model_condition.wait()
        if _initialized:
            if _load_error:
                raise Exception(_load_error)
            return _cached_model
        _loading = True
    
    try:
        if not os.path.exists(MODEL_PATH):
            raise FileNotFoundError("Model file missing at backend/models/iCW773.xml")

        # COBRApy is one of the heaviest imports in the desktop process. Keep it
        # out of the critical startup path and import it in the warm-up thread.
        import cobra

        logger.info(f"Loading SBML model from {MODEL_PATH}...")
        loaded_model = cobra.io.read_sbml_model(MODEL_PATH)
        logger.info("Model loaded successfully!")

        # ── Apply thermodynamic directionality pruning ──────────────────────
        if _THERMO_AVAILABLE:
            logger.info("Applying thermodynamic directionality pruning...")
            loaded_model, report = apply_thermodynamic_pruning(loaded_model)
            logger.info(
                f"Pruning: {report.get('n_pruned', 0)} reactions direction-locked "
                f"(forward={report.get('n_forward_locked', 0)}, "
                f"reverse={report.get('n_reverse_locked', 0)}), "
                f"data coverage={report.get('data_coverage_pct', 0)}%"
            )
        else:
            logger.info("Thermodynamic pruning skipped (thermo_pruner not available).")

        with _model_condition:
            _cached_model = loaded_model
            _load_error = None
            _initialized = True
        return loaded_model
    except Exception as e:
        error = str(e) if isinstance(e, FileNotFoundError) else f"Failed to parse SBML model: {str(e)}"
        with _model_condition:
            _load_error = error
            _initialized = True
        logger.error(error)
        raise Exception(error)
    finally:
        with _model_condition:
            _loading = False
            _model_condition.notify_all()


def start_model_preload(delay_seconds=0.01):
    """Warm the model cache without delaying API and desktop UI readiness."""
    global _preload_thread
    with _model_condition:
        if _initialized or _loading or (_preload_thread and _preload_thread.is_alive()):
            return _preload_thread
        _preload_thread = threading.Timer(
            delay_seconds,
            _preload_model_safely,
        )
        _preload_thread.name = "metabolic-model-preload"
        _preload_thread.daemon = True
        _preload_thread.start()
        return _preload_thread


def _preload_model_safely():
    try:
        load_model_if_needed()
    except Exception as exc:
        logger.warning("Background model preload failed; requests will report the cached error: %s", exc)
