import json
from eos_workflow.inspect import eos_inspect
from eos_workflow.delta_metric import birch_murnaghan_fit


def obtain_v0(filepath="eos_fitting_results.json"):
    with open(filepath, 'r') as fp:
        eos_results = json.load(fp)
    for element, results in eos_results.items():
        for config, res in results.items():
            result = birch_murnaghan_fit(res)
            energy0 = result["energy0"]
            print(energy0)


if __name__ == "__main__":
    # obtain_v0()
    eos_inspect()

