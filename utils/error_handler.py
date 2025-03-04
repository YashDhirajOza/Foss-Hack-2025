class ChemicalError(Exception):
    """Base class for chemical-related errors"""
    pass

class ReactionError(ChemicalError):
    """Error in reaction prediction"""
    pass

class ValidationError(ChemicalError):
    """Error in input validation"""
    pass

def handle_chemical_error(error):
    """Handle chemical-related errors"""
    error_types = {
        ReactionError: "Reaction calculation failed",
        ValidationError: "Invalid input parameters",
        ChemicalError: "General chemical error"
    }
    
    error_type = type(error)
    base_message = error_types.get(error_type, "Unknown error")
    return {
        "success": False,
        "error_type": error_type.__name__,
        "message": f"{base_message}: {str(error)}"
    }
