import pytest
from utils.mechanism_viewer import draw_reaction_mechanism
from data.reaction_predictor import ReactionPredictor

def test_acid_base_reaction():
    predictor = ReactionPredictor()
    result = predictor.predict_product(
        {"formula": "HCl"},
        {"formula": "NaOH"},
        {"temperature": 25, "pressure": 1}
    )
    assert result["success"] == True
    assert "NaCl" in result["products"]
    assert "H2O" in result["products"]

def test_invalid_reaction():
    predictor = ReactionPredictor()
    result = predictor.predict_product(
        {"formula": "InvalidCompound"},
        {"formula": "NaOH"},
        {"temperature": 25, "pressure": 1}
    )
    assert result["success"] == False
