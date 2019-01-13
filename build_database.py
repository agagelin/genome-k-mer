import json
import os
from pathlib import Path

import progressbar

data_dir = Path("./data/running_window_data")
files = os.listdir(str(data_dir))
output_file = Path("./data/running_window_fulldataset.json")

data = dict()
for json_file in files:
    with open(str(data_dir / json_file), "r") as f:
        current = json.load(f)
    name = current["name"]
    if name not in data:
        data[name] = current
    else:
        for k in current["distances"]:
            if k not in data[name]["distances"]:
                data[name]["distances"][k] = current["distances"][k]
            else:
                data[name]["distance"][k].update(current["distances"][k])
data = [data[k] for k in data]

with open(str(output_file), "w") as f:
    json.dump(data, f)
