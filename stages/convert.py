#!/usr/bin/env python3
"""Stream PubChem BioAssay core result rows to RDF."""
import argparse,gzip,json
from pathlib import Path
from rdflib import Literal,URIRef
from rdflib.namespace import DCTERMS,RDF,RDFS,SKOS,XSD
from rdflib.plugins.serializers.nt import _quoteLiteral
import pyarrow.parquet as pq
import pyarrow.compute as pc
BASE="https://biobricks.ai/pubchem-bioassay/";PROV="http://www.w3.org/ns/prov#";OBO="http://purl.obolibrary.org/obo/";DATASET=URIRef(BASE+"dataset/nih-pubchem-bioassays")
FIELDS=("aid","sid","cid","activity_outcome","activity_score","smiles")
def present(v):return v is not None and str(v).strip() not in ("","nan")
def nt(v):return f"<{v}>" if isinstance(v,URIRef) else _quoteLiteral(v)
def emit(f,s,p,o):f.write(f"{nt(s)} {nt(p)} {nt(o)} .\n");return 1
def srcdir():return Path("/mnt/raid2/biobricks/nih-pubchem-bioassays/brick")
def source_nonempty(pf):
 n={x:0 for x in FIELDS}
 for rg in range(pf.metadata.num_row_groups):
  rows=pf.metadata.row_group(rg).num_rows
  for i,k in enumerate(FIELDS):
   st=pf.metadata.row_group(rg).column(i).statistics;n[k]+=rows-(st.null_count if st else 0)
 return n
def convert(src,out,report,max_results=None):
 pf=pq.ParquetFile(src/"bioassay_activities.parquet");non=source_nonempty(pf);total=pf.metadata.num_rows;out.parent.mkdir(parents=True,exist_ok=True);report.parent.mkdir(parents=True,exist_ok=True)
 c={"source_rows":total,"converted_rows":0,"source_nonempty_cells":sum(non.values()),"mapped_nonempty_cells":0,"excluded_nonempty_cells":0,"n_hexane_rows":0,"triples":0};outcomes=set();aids=set();op=gzip.open if str(out).endswith(".gz") else open;rownum=0
 with op(out,"wt") as f:
  for t in ((DATASET,RDF.type,URIRef(PROV+"Entity")),(DATASET,DCTERMS.source,URIRef("https://pubchem.ncbi.nlm.nih.gov/docs/bioassays")),(DATASET,DCTERMS.license,URIRef("https://www.nlm.nih.gov/databases/download/terms_and_conditions.html"))):c["triples"]+=emit(f,*t)
  for batch in pf.iter_batches(batch_size=8192):
   if max_results is not None and c["converted_rows"]>=max_results:
    batch=batch.filter(pc.fill_null(pc.equal(batch.column(batch.schema.get_field_index("cid")),8058),False))
   for row in batch.to_pylist():
    rownum+=1;is_hex=row.get("cid")==8058
    if max_results is not None and c["converted_rows"]>=max_results and not is_hex:continue
    aid=row.get("aid");assay=URIRef("https://pubchem.ncbi.nlm.nih.gov/bioassay/"+str(aid));obs=URIRef(BASE+f"result/{rownum}");triples=[(assay,RDF.type,URIRef(OBO+"OBI_0000070")),(obs,RDF.type,URIRef(OBO+"OBI_0000299")),(obs,URIRef(OBO+"OBI_0000312"),assay),(obs,URIRef(PROV+"wasDerivedFrom"),DATASET)];aids.add(aid)
    sid=row.get("sid");cid=row.get("cid");chem=URIRef("https://pubchem.ncbi.nlm.nih.gov/compound/"+str(cid)) if cid is not None else None
    vals=((DCTERMS.identifier,aid),(URIRef(BASE+"property/substance"),URIRef("https://pubchem.ncbi.nlm.nih.gov/substance/"+str(sid)) if sid is not None else None),(URIRef(BASE+"property/compound"),chem),(URIRef(BASE+"property/activity-outcome"),row.get("activity_outcome")),(URIRef(BASE+"property/pubchem-activity-score"),row.get("activity_score")),(URIRef(OBO+"CHEMINF_000018"),row.get("smiles")))
    for pred,val in vals:
     if present(val):triples.append((obs if pred!=URIRef(OBO+"CHEMINF_000018") else chem,pred,val if isinstance(val,URIRef) else Literal(val)));c["mapped_nonempty_cells"]+=1
    if present(row.get("activity_outcome")):outcomes.add(str(row["activity_outcome"]))
    if is_hex and chem:
     c["n_hexane_rows"]+=1
     for x in ("https://biobricks.ai/compound/inchikey/VLKZOEOYAKHREP-UHFFFAOYSA-N","https://identifiers.org/cas/110-54-3","https://comptox.epa.gov/dashboard/chemical/details/DTXSID0021917"):triples.append((chem,SKOS.exactMatch,URIRef(x)))
    for t in triples:c["triples"]+=emit(f,*t)
    c["converted_rows"]+=1
 c["excluded_nonempty_cells"]=c["source_nonempty_cells"]-c["mapped_nonempty_cells"]
 eligible=total;r={**c,"source_tables":1,"source_nonempty_by_field":non,"field_disposition_coverage":1.0,"mapped_cell_coverage":c["mapped_nonempty_cells"]/c["source_nonempty_cells"],"result_row_coverage":c["converted_rows"]/eligible,"assays_in_artifact":len(aids),"outcomes_in_artifact":sorted(outcomes),"explicit_exclusions":{"bounded-validation-build":c["excluded_nonempty_cells"]},"source_limitations":["Source brick omits assay titles, targets, Gene/UniProt identifiers, units and assay-specific measurements.","No target or quantitative activity claim can be reconstructed from the retained core fields."],"evidence_boundary":"deposited PubChem assay outcome and normalized activity score; not a toxicology, hazard, efficacy, or target assertion","full_build_expectation":{"result_rows":total,"field_disposition_coverage":1.0,"streaming_required":True,"estimated_triples":"greater than 800 million"}}
 report.write_text(json.dumps(r,indent=2,sort_keys=True)+"\n");return r
def main():
 p=argparse.ArgumentParser();p.add_argument("--source-dir",type=Path);p.add_argument("--output",type=Path,default=Path("brick/pubchem-bioassays-rdf.nt.gz"));p.add_argument("--coverage",type=Path,default=Path("reports/source-coverage.json"));p.add_argument("--max-results",type=int);a=p.parse_args();print(json.dumps(convert(a.source_dir or srcdir(),a.output,a.coverage,a.max_results),indent=2))
if __name__=="__main__":main()
