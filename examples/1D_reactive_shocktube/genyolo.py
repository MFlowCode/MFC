import mfc.viz
from tqdm import tqdm

case = mfc.viz.Case(".")
lines = []

for var in tqdm(range(1, 15+1), desc="Loading and Generating Arrays"):
    case.load_variable(f"{var}", f"prim.{var}")
    elems = [str(x) for x in case.get_data()[59760][str(var)]]
    lines.append(f"real(kind(0d0)) :: var{var}(0:400) = [ &")
    for i, element in enumerate(elems):
        if i == len(elems) - 1:
            lines.append(f"{element} &")
        else:
            lines.append(f"{element}, &")
    lines.append("]")

for var in tqdm(range(1, 15+1), desc="Loading Conserved Variables"):
    lines.append(f"q_prim_vf({var})%sf(i, j, 0) = var{var}(i)")

print("\n".join(lines))