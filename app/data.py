from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from html import unescape
from typing import Any

import httpx


HTTP_TIMEOUT = 12.0
USER_AGENT = "BioLens/0.1 research-prototype"


@dataclass
class ParsedQuery:
    raw: str
    gene: str
    mutation: str


def parse_query(query: str) -> ParsedQuery:
    parts = query.replace(",", " ").split()
    gene = parts[0].upper() if parts else ""
    mutation = " ".join(parts[1:]) if len(parts) > 1 else ""
    return ParsedQuery(raw=query, gene=gene, mutation=mutation)


def build_analysis(query: str) -> dict[str, Any]:
    parsed = parse_query(query)
    with httpx.Client(timeout=HTTP_TIMEOUT, headers={"User-Agent": USER_AGENT}) as client:
        mygene = safe_fetch(fetch_mygene, client, parsed.gene, fallback={})
        uniprot = safe_fetch(fetch_uniprot, client, parsed.gene, fallback={})
        clinvar = safe_fetch(fetch_clinvar, client, parsed, fallback={"diseases": []})
        papers = safe_fetch(fetch_pubmed, client, parsed, fallback=[])
        pathways = safe_fetch(fetch_reactome_pathways, client, uniprot.get("accession"), fallback=[])

    gene_symbol = mygene.get("symbol") or parsed.gene
    mutation = parsed.mutation or "No mutation specified"
    accession = uniprot.get("accession")
    clinvar_location = clinvar.get("location", {})

    biological_processes = unique(
        uniprot.get("biological_processes", []) + mygene.get("go_bp", [])
    )[:6]
    molecular_functions = unique(
        uniprot.get("molecular_functions", []) + mygene.get("go_mf", [])
    )[:6]

    overview_refs = []
    if accession:
        overview_refs.append(f"UniProt:{accession}")
    if clinvar.get("accession"):
        overview_refs.append(f"ClinVar:{clinvar['accession']}")
    if mygene.get("entrezgene"):
        overview_refs.append(f"NCBI Gene:{mygene['entrezgene']}")
    if not overview_refs:
        overview_refs.append("No source identifiers returned")

    return {
        "query": query,
        "overview": {
            "gene": gene_symbol or "Unknown gene",
            "mutation": mutation,
            "proteinPosition": clinvar.get("protein_change") or infer_protein_change(mutation),
            "mutationType": clinvar.get("variant_type") or "Not reported by source",
            "clinicalSignificance": clinvar.get("clinical_significance") or "Not reported by source",
            "chromosome": clinvar_location.get("band") or mygene.get("chromosome") or "Not reported by source",
            "coordinates": clinvar_location.get("coordinates") or mygene.get("coordinates") or "Not reported by source",
            "referenceIds": overview_refs,
        },
        "gene": {
            "description": mygene.get("summary") or "No gene summary returned by MyGene.info.",
            "aliases": mygene.get("aliases") or [],
            "chromosome": mygene.get("chromosome") or "Not reported by source",
            "proteinCoding": mygene.get("type_of_gene") == "protein-coding",
            "expression": "Expression data not requested from a live source in this build.",
            "length": mygene.get("length") or 0,
            "functions": biological_processes[:5] or ["No GO biological-process terms returned by source."],
        },
        "protein": {
            "name": uniprot.get("protein_name") or "No reviewed UniProt human protein returned.",
            "length": uniprot.get("length") or 0,
            "domains": uniprot.get("domains") or [],
            "subcellularLocation": uniprot.get("subcellular_location") or "Not reported by source",
            "biologicalProcesses": biological_processes or ["Not reported by source"],
            "molecularFunctions": molecular_functions or ["Not reported by source"],
            "family": uniprot.get("family") or "Not reported by source",
            "structureUrl": alphafold_url(accession),
            "accession": accession or "",
        },
        "diseases": clinvar.get("diseases") or [],
        "papers": papers,
        "aiSummary": build_deferred_summary(gene_symbol, mutation, overview_refs),
        "pathway": build_pathway_nodes(gene_symbol, uniprot, pathways, clinvar),
    }


def build_suggestions(query: str) -> list[dict[str, str]]:
    parsed = parse_query(query)
    if not parsed.gene:
        return []
    with httpx.Client(timeout=HTTP_TIMEOUT, headers={"User-Agent": USER_AGENT}) as client:
        if parsed.mutation:
            return safe_fetch(fetch_variant_suggestions, client, parsed, fallback=[])
        return safe_fetch(fetch_gene_suggestions, client, parsed.gene, fallback=[])


def fetch_gene_suggestions(client: httpx.Client, term: str) -> list[dict[str, str]]:
    response = client.get(
        "https://mygene.info/v3/query",
        params={
            "q": f"symbol:{term}* OR alias:{term}*",
            "species": "human",
            "fields": "symbol,name,alias,type_of_gene,entrezgene",
            "size": 8,
        },
    )
    response.raise_for_status()
    suggestions = []
    seen = set()
    for hit in response.json().get("hits", []):
        symbol = hit.get("symbol")
        if not symbol or symbol in seen:
            continue
        seen.add(symbol)
        gene_type = hit.get("type_of_gene") or "gene"
        suggestions.append({
            "value": symbol,
            "label": symbol,
            "detail": hit.get("name") or gene_type,
            "source": f"MyGene.info / NCBI Gene {hit.get('entrezgene')}" if hit.get("entrezgene") else "MyGene.info",
            "kind": "gene",
        })
    return suggestions


def fetch_variant_suggestions(client: httpx.Client, parsed: ParsedQuery) -> list[dict[str, str]]:
    term = " ".join(part for part in [parsed.gene, parsed.mutation] if part)
    search = client.get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
        params={
            "db": "clinvar",
            "term": term,
            "retmode": "json",
            "retmax": 12,
        },
    )
    search.raise_for_status()
    ids = search.json().get("esearchresult", {}).get("idlist", [])
    if not ids:
        protein = one_letter_to_three_letter(parsed.mutation)
        if protein:
            search = client.get(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                params={
                    "db": "clinvar",
                    "term": f"{parsed.gene} {protein}",
                    "retmode": "json",
                    "retmax": 12,
                },
            )
            search.raise_for_status()
            ids = search.json().get("esearchresult", {}).get("idlist", [])
    if not ids:
        return []

    summary = client.get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
        params={"db": "clinvar", "id": ",".join(ids), "retmode": "json"},
    )
    summary.raise_for_status()
    result = summary.json().get("result", {})
    suggestions = []
    seen = set()
    for uid in ids:
        record = result.get(uid, {})
        if not record:
            continue
        title = clean_text(record.get("title") or "")
        variation = (record.get("variation_set") or [{}])[0]
        variant = display_variant_from_record(record, parsed.mutation)
        if not variant:
            continue
        value = f"{parsed.gene} {variant}"
        key = normalize_variant_text(value)
        if key in seen:
            continue
        seen.add(key)
        classification = record.get("germline_classification") or {}
        significance = classification.get("description") or "clinical significance not reported"
        suggestions.append({
            "value": value,
            "label": value,
            "detail": title or variation.get("variation_name") or significance,
            "source": f"ClinVar {record.get('accession')}" if record.get("accession") else "ClinVar",
            "kind": "variant",
        })
        if len(suggestions) == 8:
            break
    return suggestions


def fetch_mygene(client: httpx.Client, gene: str) -> dict[str, Any]:
    if not gene:
        return {}
    response = client.get(
        "https://mygene.info/v3/query",
        params={
            "q": f"symbol:{gene}",
            "species": "human",
            "fields": "symbol,name,summary,alias,genomic_pos,type_of_gene,go,uniprot,entrezgene",
            "size": 1,
        },
    )
    response.raise_for_status()
    hits = response.json().get("hits", [])
    if not hits:
        return {}
    hit = hits[0]
    genomic = hit.get("genomic_pos") or {}
    start = genomic.get("start")
    end = genomic.get("end")
    chromosome = str(genomic.get("chr", "")) if genomic.get("chr") else ""
    aliases = hit.get("alias") or []
    if isinstance(aliases, str):
        aliases = [aliases]

    return {
        "symbol": hit.get("symbol"),
        "summary": hit.get("summary"),
        "aliases": aliases[:12],
        "chromosome": f"chr{chromosome}" if chromosome else "",
        "coordinates": f"chr{chromosome}:{start}-{end}" if chromosome and start and end else "",
        "length": abs(int(end) - int(start)) + 1 if start and end else 0,
        "type_of_gene": hit.get("type_of_gene"),
        "entrezgene": str(hit.get("entrezgene", "")),
        "go_bp": go_terms(hit, "BP"),
        "go_mf": go_terms(hit, "MF"),
    }


def fetch_uniprot(client: httpx.Client, gene: str) -> dict[str, Any]:
    if not gene:
        return {}
    response = client.get(
        "https://rest.uniprot.org/uniprotkb/search",
        params={
            "query": f"(gene_exact:{gene}) AND (organism_id:9606) AND (reviewed:true)",
            "format": "json",
            "size": 1,
        },
    )
    response.raise_for_status()
    results = response.json().get("results", [])
    if not results:
        return {}
    entry = results[0]
    comments = entry.get("comments", [])
    features = entry.get("features", [])

    return {
        "accession": entry.get("primaryAccession"),
        "protein_name": protein_name(entry),
        "length": (entry.get("sequence") or {}).get("length", 0),
        "domains": domains_from_features(features),
        "subcellular_location": subcellular_location(comments),
        "biological_processes": keywords(entry, "Biological process"),
        "molecular_functions": keywords(entry, "Molecular function"),
        "family": protein_family(comments),
    }


def fetch_clinvar(client: httpx.Client, parsed: ParsedQuery) -> dict[str, Any]:
    if not parsed.gene or not parsed.mutation:
        return {"diseases": []}
    first = {}
    for term in clinvar_search_terms(parsed):
        search = client.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
            params={
                "db": "clinvar",
                "term": term,
                "retmode": "json",
                "retmax": 10,
            },
        )
        search.raise_for_status()
        ids = search.json().get("esearchresult", {}).get("idlist", [])
        if not ids:
            continue

        summary = client.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
            params={"db": "clinvar", "id": ",".join(ids), "retmode": "json"},
        )
        summary.raise_for_status()
        result = summary.json().get("result", {})
        first = best_clinvar_record(result, ids, parsed)
        if first:
            break
    if not first:
        return {"diseases": []}
    classification = first.get("germline_classification") or {}
    variation = (first.get("variation_set") or [{}])[0]

    return {
        "accession": first.get("accession"),
        "clinical_significance": classification.get("description"),
        "review_status": classification.get("review_status"),
        "variant_type": variation.get("variant_type") or first.get("obj_type"),
        "protein_change": protein_change_from_title(first.get("title") or variation.get("variation_name", "")),
        "location": clinvar_location(variation),
        "diseases": diseases_from_clinvar(first, parsed.mutation),
    }


def fetch_pubmed(client: httpx.Client, parsed: ParsedQuery) -> list[dict[str, Any]]:
    if not parsed.gene:
        return []
    term = f"{parsed.gene} {parsed.mutation}".strip()
    response = client.get(
        "https://www.ebi.ac.uk/europepmc/webservices/rest/search",
        params={
            "query": f'({term}) AND SRC:MED',
            "format": "json",
            "pageSize": 6,
            "resultType": "core",
        },
    )
    response.raise_for_status()
    results = response.json().get("resultList", {}).get("result", [])
    papers = []
    for item in results:
        pmid = item.get("pmid") or item.get("id")
        title = clean_text(item.get("title") or "Untitled PubMed record")
        abstract = clean_text(item.get("abstractText") or "No abstract returned by Europe PMC.")
        papers.append({
            "title": title,
            "authors": split_authors(item.get("authorString")),
            "journal": item.get("journalTitle") or item.get("journalInfo", {}).get("journal", {}).get("title") or "Journal not reported",
            "year": int(item.get("pubYear") or 0),
            "abstract": abstract,
            "relevance": pubmed_relevance(title, abstract, term),
            "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else item.get("fullTextUrlList", {}).get("fullTextUrl", [{}])[0].get("url", "https://pubmed.ncbi.nlm.nih.gov/"),
        })
    return papers


def fetch_reactome_pathways(client: httpx.Client, accession: str | None) -> list[str]:
    if not accession:
        return []
    response = client.get(
        f"https://reactome.org/ContentService/data/mapping/UniProt/{accession}/pathways",
        params={"species": 9606},
    )
    if response.status_code == 404:
        return []
    response.raise_for_status()
    pathways = response.json()
    return unique([item.get("displayName", "") for item in pathways if item.get("displayName")])[:4]


def go_terms(hit: dict[str, Any], category: str) -> list[str]:
    go = (hit.get("go") or {}).get(category, [])
    if isinstance(go, dict):
        go = [go]
    return unique([item.get("term", "") for item in go if item.get("term")])[:8]


def protein_name(entry: dict[str, Any]) -> str:
    description = entry.get("proteinDescription") or {}
    recommended = description.get("recommendedName") or {}
    full = recommended.get("fullName") or {}
    return full.get("value", "")


def domains_from_features(features: list[dict[str, Any]]) -> list[dict[str, Any]]:
    domains = []
    accepted = {"Domain", "Region", "Repeat", "Motif", "Topological domain"}
    for feature in features:
        if feature.get("type") not in accepted:
            continue
        location = feature.get("location") or {}
        start = (((location.get("start") or {}).get("value")) or 0)
        end = (((location.get("end") or {}).get("value")) or start)
        description = feature.get("description") or feature.get("type")
        if start and end and description:
            domains.append({"name": description[:60], "start": int(start), "end": int(end)})
        if len(domains) == 6:
            break
    return domains


def subcellular_location(comments: list[dict[str, Any]]) -> str:
    locations = []
    for comment in comments:
        if comment.get("commentType") != "SUBCELLULAR LOCATION":
            continue
        for item in comment.get("subcellularLocations", []):
            location = (item.get("location") or {}).get("value")
            if location:
                locations.append(location)
    return "; ".join(unique(locations)[:5])


def protein_family(comments: list[dict[str, Any]]) -> str:
    for comment in comments:
        if comment.get("commentType") != "SIMILARITY":
            continue
        texts = comment.get("texts") or []
        if texts and texts[0].get("value"):
            return texts[0]["value"]
    return ""


def keywords(entry: dict[str, Any], category: str) -> list[str]:
    return unique([
        keyword.get("name", "")
        for keyword in entry.get("keywords", [])
        if keyword.get("category") == category and keyword.get("name")
    ])[:8]


def clinvar_location(variation: dict[str, Any]) -> dict[str, str]:
    current = {}
    for loc in variation.get("variation_loc", []):
        if loc.get("assembly_name") == "GRCh38":
            current = loc
            break
    if not current and variation.get("variation_loc"):
        current = variation["variation_loc"][0]
    if not current:
        return {}
    chr_name = current.get("chr")
    start = current.get("display_start") or current.get("start")
    stop = current.get("display_stop") or current.get("stop")
    assembly = current.get("assembly_name") or "assembly"
    coordinates = f"{assembly} chr{chr_name}:{start}-{stop}" if chr_name and start and stop else ""
    return {"band": current.get("band", ""), "coordinates": coordinates}


def diseases_from_clinvar(record: dict[str, Any], mutation: str) -> list[dict[str, Any]]:
    trait_set = record.get("trait_set") or []
    classification = record.get("germline_classification") or {}
    review = classification.get("review_status") or "ClinVar review status not reported"
    significance = classification.get("description") or "Clinical significance not reported"
    diseases = []
    for index, trait in enumerate(trait_set[:6]):
        name = trait.get("trait_name") or trait.get("name")
        if not name:
            continue
        diseases.append({
            "name": name,
            "description": f"ClinVar trait associated with {record.get('title', 'this variant')}.",
            "evidence": f"{significance}; {review}.",
            "mutation": mutation or record.get("title", ""),
            "confidence": max(0.55, 0.95 - index * 0.06),
            "severity": severity_for_significance(significance),
            "link": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{record.get('uid', '')}/",
        })
    if diseases:
        return diseases
    if record.get("title"):
        return [{
            "name": record["title"],
            "description": "ClinVar returned a matching variant but no trait list in the summary response.",
            "evidence": f"{significance}; {review}.",
            "mutation": mutation or record["title"],
            "confidence": 0.72,
            "severity": severity_for_significance(significance),
            "link": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{record.get('uid', '')}/",
        }]
    return []


def display_variant_from_record(record: dict[str, Any], typed_mutation: str) -> str:
    title = record.get("title") or ""
    variation = (record.get("variation_set") or [{}])[0]
    cdna = variation.get("cdna_change") or ""
    protein = protein_change_from_title(title)

    typed = typed_mutation.strip()
    if typed.lower().startswith("c.") and cdna:
        return cdna
    if typed.lower().startswith("p.") and protein:
        return protein
    if protein:
        return protein
    if cdna:
        return cdna

    aliases = variation.get("aliases") or []
    if aliases:
        return aliases[0]
    name = variation.get("variation_name") or title
    match = re.search(r":(c\.[^( ]+)", name)
    return match.group(1) if match else ""


def best_clinvar_record(result: dict[str, Any], ids: list[str], parsed: ParsedQuery) -> dict[str, Any]:
    records = [result.get(uid, {}) for uid in ids if result.get(uid)]
    if not parsed.mutation:
        return records[0] if records else {}
    for record in records:
        if clinvar_record_matches(record, parsed):
            return record
    return {}


def clinvar_search_terms(parsed: ParsedQuery) -> list[str]:
    terms = [f"{parsed.gene} {parsed.mutation}"]
    protein = one_letter_to_three_letter(parsed.mutation)
    if protein:
        terms.append(f"{parsed.gene} {protein}")
        terms.append(f"{parsed.gene} {protein.removeprefix('p.')}")
    return unique(terms)


def clinvar_record_matches(record: dict[str, Any], parsed: ParsedQuery) -> bool:
    mutation_tokens = variant_tokens(parsed.mutation)
    if not mutation_tokens:
        return False
    searchable = " ".join([
        record.get("title", ""),
        record.get("accession", ""),
        " ".join((record.get("variation_set") or [{}])[0].get("aliases") or []),
        (record.get("variation_set") or [{}])[0].get("variation_name", ""),
        (record.get("variation_set") or [{}])[0].get("cdna_change", ""),
    ])
    normalized = normalize_variant_text(searchable)
    return any(token in normalized for token in mutation_tokens)


def variant_tokens(mutation: str) -> list[str]:
    normalized = normalize_variant_text(mutation)
    tokens = [normalized]
    protein = one_letter_to_three_letter(mutation)
    if protein:
        tokens.append(normalize_variant_text(protein))
    if mutation.lower().startswith("c."):
        tokens.append(normalize_variant_text(mutation[2:]))
    return unique(tokens)


AA3 = {
    "A": "Ala",
    "R": "Arg",
    "N": "Asn",
    "D": "Asp",
    "C": "Cys",
    "Q": "Gln",
    "E": "Glu",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "L": "Leu",
    "K": "Lys",
    "M": "Met",
    "F": "Phe",
    "P": "Pro",
    "S": "Ser",
    "T": "Thr",
    "W": "Trp",
    "Y": "Tyr",
    "V": "Val",
}


def one_letter_to_three_letter(mutation: str) -> str:
    match = re.fullmatch(r"p?\.?([A-Z])(\d+)([A-Z])", mutation.strip())
    if not match:
        return ""
    ref, pos, alt = match.groups()
    if ref not in AA3 or alt not in AA3:
        return ""
    return f"p.{AA3[ref]}{pos}{AA3[alt]}"


def normalize_variant_text(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.lower())


def parse_pubmed_articles(xml: str, term: str) -> list[dict[str, Any]]:
    root = ET.fromstring(xml)
    papers = []
    for article in root.findall(".//PubmedArticle"):
        medline = article.find("./MedlineCitation")
        pmid = text(medline, "./PMID") if medline is not None else ""
        article_node = medline.find("./Article") if medline is not None else None
        if article_node is None:
            continue
        title = "".join(article_node.findtext("./ArticleTitle", default="").split())
        title = article_node.findtext("./ArticleTitle", default="Untitled PubMed record")
        journal = article_node.findtext("./Journal/Title", default="Journal not reported")
        year = first_year(article_node)
        authors = []
        for author in article_node.findall("./AuthorList/Author")[:5]:
            last = author.findtext("./LastName", default="")
            initials = author.findtext("./Initials", default="")
            collective = author.findtext("./CollectiveName", default="")
            name = collective or " ".join(part for part in [initials, last] if part)
            if name:
                authors.append(name)
        abstract_parts = [
            "".join(part.itertext())
            for part in article_node.findall("./Abstract/AbstractText")
        ]
        papers.append({
            "title": title,
            "authors": authors or ["Authors not reported"],
            "journal": journal,
            "year": year,
            "abstract": " ".join(abstract_parts) or "No abstract returned by PubMed.",
            "relevance": pubmed_relevance(title, " ".join(abstract_parts), term),
            "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else "https://pubmed.ncbi.nlm.nih.gov/",
        })
    return papers


def split_authors(author_string: str | None) -> list[str]:
    if not author_string:
        return ["Authors not reported"]
    authors = [author.strip().rstrip(".") for author in author_string.split(",") if author.strip()]
    return authors[:6] or ["Authors not reported"]


def clean_text(value: str) -> str:
    without_tags = re.sub(r"<[^>]+>", "", unescape(value))
    return " ".join(without_tags.split())


def build_deferred_summary(gene: str, mutation: str, refs: list[str]) -> dict[str, Any]:
    return {
        "summary": "AI synthesis is intentionally deferred. The surrounding dashboard sections are populated from live public data sources.",
        "mutationImpact": f"No model-generated interpretation is being produced for {gene} {mutation}. Review the retrieved source records directly.",
        "biologicalEffects": [
            "Use the live UniProt, ClinVar, Reactome, and PubMed evidence shown elsewhere on this page.",
            "No unstated biological conclusions are generated in this placeholder.",
        ],
        "notes": ["Generated summary disabled by request.", "Facts shown outside this box come from live public APIs when available."],
        "confidence": 0,
        "citations": refs,
    }


def build_pathway_nodes(gene: str, uniprot: dict[str, Any], pathways: list[str], clinvar: dict[str, Any]) -> list[dict[str, str]]:
    protein_label = uniprot.get("protein_name") or uniprot.get("accession") or "Protein not returned"
    pathway_label = pathways[0] if pathways else "Reactome pathway not returned"
    disease_label = (clinvar.get("diseases") or [{}])[0].get("name", "ClinVar trait not returned")
    return [
        {"id": "gene", "label": gene or "Gene", "kind": "gene"},
        {"id": "protein", "label": protein_label[:36], "kind": "protein"},
        {"id": "pathway", "label": pathway_label[:42], "kind": "pathway"},
        {"id": "disease", "label": disease_label[:42], "kind": "disease"},
        {"id": "source", "label": "Live public evidence", "kind": "source"},
    ]


def alphafold_url(accession: str | None) -> str:
    if not accession:
        return ""
    return f"https://molstar.org/viewer/?url=https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v4.cif&type=mmcif"


def infer_protein_change(mutation: str) -> str:
    if mutation.startswith("p."):
        return mutation
    match = re.fullmatch(r"([A-Z])(\d+)([A-Z]|del|fs.*)?", mutation)
    if match:
        return f"p.{mutation}"
    return "Not reported by source"


def protein_change_from_title(title: str) -> str:
    match = re.search(r"\((p\.[^)]+)\)", title)
    return match.group(1) if match else ""


def severity_for_significance(significance: str) -> str:
    lowered = significance.lower()
    if "pathogenic" in lowered:
        return "critical"
    if "risk" in lowered or "drug" in lowered or "conflicting" in lowered:
        return "high"
    return "moderate"


def pubmed_relevance(title: str, abstract: str, term: str) -> float:
    haystack = f"{title} {abstract}".lower()
    tokens = [token.lower() for token in term.split() if token]
    if not tokens:
        return 0.75
    hits = sum(1 for token in tokens if token in haystack)
    return min(0.99, 0.65 + hits * 0.12)


def first_year(article_node: ET.Element) -> int:
    for path in [
        "./Journal/JournalIssue/PubDate/Year",
        "./ArticleDate/Year",
        "./Journal/JournalIssue/PubDate/MedlineDate",
    ]:
        value = article_node.findtext(path, default="")
        match = re.search(r"\d{4}", value)
        if match:
            return int(match.group(0))
    return 0


def text(node: ET.Element | None, path: str) -> str:
    if node is None:
        return ""
    return node.findtext(path, default="")


def unique(items: list[str]) -> list[str]:
    seen = set()
    values = []
    for item in items:
        clean = " ".join(str(item).split())
        if clean and clean.lower() not in seen:
            seen.add(clean.lower())
            values.append(clean)
    return values


def safe_fetch(function, *args, fallback):
    try:
        return function(*args)
    except (httpx.HTTPError, ValueError, ET.ParseError):
        return fallback
