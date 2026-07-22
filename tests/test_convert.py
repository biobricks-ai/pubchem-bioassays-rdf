import gzip,importlib.util
from pathlib import Path
ROOT=Path(__file__).parents[1];s=importlib.util.spec_from_file_location("c",ROOT/"stages/convert.py");M=importlib.util.module_from_spec(s);s.loader.exec_module(M);SRC=Path("/mnt/raid2/biobricks/nih-pubchem-bioassays/brick")
def test_streaming_real_source_and_hexane(tmp_path):
 out=tmp_path/"x.nt.gz";r=M.convert(SRC,out,tmp_path/"r.json",100)
 assert r["source_rows"]==92224748 and r["converted_rows"]==105 and r["n_hexane_rows"]==5
 assert r["field_disposition_coverage"]==1 and 0<r["result_row_coverage"]<1
 with gzip.open(out,"rt") as f:t=f.read()
 assert "pubchem.ncbi.nlm.nih.gov/bioassay/1" in t and "pubchem.ncbi.nlm.nih.gov/compound/8058" in t
 assert "VLKZOEOYAKHREP-UHFFFAOYSA-N" in t and "target" not in t.lower()
def test_source_limitations_are_explicit(tmp_path):
 r=M.convert(SRC,tmp_path/"x.nt.gz",tmp_path/"r.json",1)
 assert any("Gene/UniProt" in x for x in r["source_limitations"])
 assert "not a toxicology" in r["evidence_boundary"]
