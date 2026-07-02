from fastapi import FastAPI, Query
from fastapi.middleware.cors import CORSMiddleware

from .data import build_analysis, build_suggestions
from .models import Analysis, Suggestion

app = FastAPI(title="BioLens API", version="0.1.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/health")
def health() -> dict[str, str]:
    return {"status": "ok"}


@app.get("/api/analyze", response_model=Analysis)
def analyze(q: str = Query(..., min_length=2, description="Gene or gene plus mutation")) -> Analysis:
    return Analysis.model_validate(build_analysis(q))


@app.get("/api/suggest", response_model=list[Suggestion])
def suggest(q: str = Query(..., min_length=1, description="Partial gene or gene plus mutation")) -> list[Suggestion]:
    return [Suggestion.model_validate(item) for item in build_suggestions(q)]
