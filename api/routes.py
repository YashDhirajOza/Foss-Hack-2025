from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from data.reaction_predictor import ReactionPredictor

app = FastAPI()
predictor = ReactionPredictor()

class ReactionRequest(BaseModel):
    reactant1: str
    reactant2: str
    temperature: float
    pressure: float

@app.post("/predict_reaction/")
async def predict_reaction(request: ReactionRequest):
    try:
        result = predictor.predict_product(
            {"formula": request.reactant1},
            {"formula": request.reactant2},
            {"temperature": request.temperature, "pressure": request.pressure}
        )
        return result
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))
