import re

def remove_higher_dimensional_keys(case: dict, ndims: int) -> dict:
    assert 1 <= ndims <= 3

    rm_dims = [set(), set(['y', 'z']), set(['z']), set()][ndims]
    dim_ids = {dim: i + 1 for i, dim in enumerate(['x', 'y', 'z'])}
    dim_mnp = {'x': 'm', 'y': 'n', 'z': 'p'}

    rm_keys = set()
    for key in case.keys():
        for dim in rm_dims:
            if any([
                re.match(f'.+_{dim}', key), re.match(f'{dim}_.+', key),
                re.match(f'%{dim}', key),   f'%vel({dim_ids[dim]})' in key
            ]):
                rm_keys.add(key)
                break

    new_case = {k: v for k, v in case.items() if k not in rm_keys}

    for dim in rm_dims:
        new_case[dim_mnp[dim]] = 0

    return new_case
