from pydantic import BaseModel, Field


class Overview(BaseModel):
    gene: str
    mutation: str
    proteinPosition: str
    mutationType: str
    clinicalSignificance: str
    chromosome: str
    coordinates: str
    referenceIds: list[str]


class GeneInfo(BaseModel):
    description: str
    aliases: list[str]
    chromosome: str
    proteinCoding: bool
    expression: str
    length: int
    functions: list[str]


class Domain(BaseModel):
    name: str
    start: int
    end: int


class ProteinInfo(BaseModel):
    name: str
    length: int
    domains: list[Domain]
    subcellularLocation: str
    biologicalProcesses: list[str]
    molecularFunctions: list[str]
    family: str
    structureUrl: str = ""
    accession: str = ""


class Disease(BaseModel):
    name: str
    description: str
    evidence: str
    mutation: str
    confidence: float = Field(ge=0, le=1)
    severity: str
    link: str


class Paper(BaseModel):
    title: str
    authors: list[str]
    journal: str
    year: int
    abstract: str
    relevance: float = Field(ge=0, le=1)
    url: str


class AISummary(BaseModel):
    summary: str
    mutationImpact: str
    biologicalEffects: list[str]
    notes: list[str]
    confidence: float = Field(ge=0, le=1)
    citations: list[str]


class PathwayNode(BaseModel):
    id: str
    label: str
    kind: str


class Suggestion(BaseModel):
    value: str
    label: str
    detail: str
    source: str
    kind: str


class Analysis(BaseModel):
    query: str
    overview: Overview
    gene: GeneInfo
    protein: ProteinInfo
    diseases: list[Disease]
    papers: list[Paper]
    aiSummary: AISummary
    pathway: list[PathwayNode]
