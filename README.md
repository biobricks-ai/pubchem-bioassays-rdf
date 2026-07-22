# pubchem-bioassays-rdf

Streaming RDF transformation of 92,224,748 PubChem BioAssay core result rows. It preserves canonical PubChem AID, SID and CID links, deposited activity outcome, PubChem normalized activity score, reported SMILES, and provenance. It does not reinterpret screening results as toxicology, hazard, efficacy, or target evidence.

The upstream brick retained only six uniform fields and discarded assay metadata, targets, Gene/UniProt identifiers, units, and assay-specific measurements. This graph therefore makes no target or quantitative potency assertion. Recovering those fields requires a separately versioned source brick rebuilt from the raw PubChem assay archives or APIs.

NCBI database content produced by the US Government is generally public domain in the United States, but PubChem includes depositor content that can retain third-party rights. The graph points to the NLM download terms and retains source provenance rather than claiming a Creative Commons license.

Run `python stages/convert.py --max-results 10000` and `pytest -q`. Omit the limit for a full streaming build. The bounded artifact retains every n-hexane row beyond its cap and reports row coverage separately from 100% source-field disposition. A full graph is expected to exceed 800 million triples.
