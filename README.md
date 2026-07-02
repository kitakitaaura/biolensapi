# BioLens Backend

This backend serves the BioLens API using FastAPI and Python. It fetches gene, protein, variant, and literature evidence from live public sources and exposes analysis and autosuggest endpoints for the frontend.

## Run locally

```bash
cd backend
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
uvicorn app.main:app --reload --host 0.0.0.0 --port 8000
```

## Endpoints

- `GET /api/analyze?q=<query>` — returns structured analysis for a gene or gene plus mutation.
- `GET /api/suggest?q=<query>` — returns autosuggest candidates for gene or variant search.
- `GET /health` — health check endpoint.

## Notes

- The backend uses CORS middleware so the frontend can call it from a separate origin.
- If you deploy it separately, point the frontend to the backend service URL using `NEXT_PUBLIC_API_URL`.
