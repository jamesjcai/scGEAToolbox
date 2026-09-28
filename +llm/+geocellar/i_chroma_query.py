"""
i_chroma_query.py  —  Query GEOcellar ChromaDB Cloud and print JSON results.

Usage:
    CHROMA_API_KEY=... CHROMA_QUERY_TOPIC=... python i_chroma_query.py [n_results]
    python i_chroma_query.py <topic> <chroma_api_key> [n_results]   (older form)

The key and the topic are read from the environment: on the command line
the key is visible to every process on the machine (ps, Task Manager), and
a quote in the topic broke the shell command chroma_query.m built.

Outputs a JSON array of study objects to stdout.
Requires: chromadb, sentence-transformers
"""

import json
import os
import sys

import chromadb
from chromadb.utils import embedding_functions

TENANT     = "87f81beb-4f84-4ef7-a136-977cd157d49b"
DATABASE   = "GEOcellar_Repo"
COLLECTION = "geo"


def main():
    topic = os.environ.get("CHROMA_QUERY_TOPIC")
    chroma_api_key = os.environ.get("CHROMA_API_KEY")
    if topic and chroma_api_key:
        n_results = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    elif len(sys.argv) >= 3:
        topic = sys.argv[1]
        chroma_api_key = sys.argv[2]
        n_results = int(sys.argv[3]) if len(sys.argv) > 3 else 20
    else:
        print("Usage: CHROMA_API_KEY=... CHROMA_QUERY_TOPIC=... "
              "i_chroma_query.py [n_results]", file=sys.stderr)
        sys.exit(1)

    client = chromadb.HttpClient(
        ssl=True,
        host="api.trychroma.com",
        tenant=TENANT,
        database=DATABASE,
        headers={"x-chroma-token": chroma_api_key},
    )

    ef = embedding_functions.SentenceTransformerEmbeddingFunction(
        model_name="all-MiniLM-L6-v2"
    )
    collection = client.get_collection(COLLECTION, embedding_function=ef)

    results = collection.query(
        query_texts=[topic],
        n_results=n_results,
        include=["metadatas", "documents", "distances"],
    )

    studies = []
    metas = results.get("metadatas", [[]])[0]
    docs  = results.get("documents", [[]])[0]

    ids = results.get("ids", [[]])[0]
    for i, m in enumerate(metas):
        doc_id = ids[i] if i < len(ids) else ""
        studies.append({
            "id":             m.get("id", doc_id),
            # Sample records carry the series they belong to; studies do not
            "parent_study":   m.get("ParentStudy", ""),
            "title":          m.get("Title", ""),
            "summary":        (m.get("Summary") or (docs[i] if i < len(docs) else ""))[:500],
            # Samples before GSMList: Samples is comma-joined accessions,
            # GSMList is a Python list repr ("['GSM123']"), brackets included
            "samples":        m.get("Samples") or m.get("GSMList", ""),
            "organism":       m.get("Organism", ""),
            "overall_design": (m.get("OverallDesign") or "")[:300],
        })

    print(json.dumps(studies))


if __name__ == "__main__":
    main()
