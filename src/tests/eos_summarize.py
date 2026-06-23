import os.path
import json
from eos_workflow.utilities import ACWF_CONFIGURATIONS
import pandas as pd

eos_path = "/home/wjing/.abinit/pseudos/ONCVPSP-PBE-SR-PDv0.6/Ir/exc/test"
eos_result_basename = 'Ir-eos_fitting_results'
n = 10

nu = {}

for i in range(10, 9+n):
    result = os.path.join(eos_path, f"{eos_result_basename}_{i}.json")
    with open(result, 'r+') as fp:
        js = json.load(fp)
    tmp = dict()
    for el, res in js.items():
        for config, output in res.items():
            if config in ACWF_CONFIGURATIONS:
                tmp[config] = output["rel_errors_vec_length"]
    nu[i] = tmp

df = pd.DataFrame.from_dict(nu, orient="index")
df.to_csv(f"{eos_result_basename}_2.csv")


