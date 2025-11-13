from atomate2.abinit.sets.base import as_pseudo_table
import json

pseudos = "ONCVPSP-PBE-SR-PDv0.6:standard"
djson_path = '/home/wjing/.abinit/pseudos/ONCVPSP-PBE-SR-PDv0.6/standard.djson'
# pseudos = as_pseudo_table(pseudos)

with open(djson_path, "rt") as fh:
    djson = json.load(fh)
