import json
from pathlib import Path
R=Path(__file__).parents[1]
def test_health():
 c=json.loads((R/"health/source-coverage.json").read_text());o=json.loads((R/"health/ontology-policy.json").read_text())
 assert c["full_source_rows"]==92224748 and c["bounded_n_hexane_rows"]==5
 assert c["bounded_field_disposition_coverage"]==1 and c["full_coverage_status"]=="supported but not executed"
 assert o["new_classes"]==0 and o["target_assertions"]==o["toxicology_assertions"]==0
