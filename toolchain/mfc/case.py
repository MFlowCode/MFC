import json, math, copy, dataclasses

@dataclasses.dataclass(init=False)
class Case:
    params: dict

    def __init__(self, params: dict) -> None:
        self.params = copy.deepcopy(params)

    def get_parameters(self) -> str:
        return self.params.keys()

    def get_cell_count(self) -> int:
        return math.prod([max(1, int(self.params.get(dir, 0))) for dir in ["m", "n", "p"]])

    def has_parameter(self, key: str)-> bool:
        return key in self.get_parameters()

    def gen_json_dict_str(self) -> str:
        return json.dumps(self.params, indent=4)

    def __getitem__(self, key: str) -> str:
        return self.params[key]

    def __setitem__(self, key: str, val: str):
        self.params[key] = val
