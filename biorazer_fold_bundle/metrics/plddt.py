import pickle
from pathlib import Path
from biorazer.structure.io.protein import CIF2CIF


def report_atom_plddts(self, seed, sample):
    samples = self.samples
    selector = (samples["seed"] == seed) & (samples["sample"] == sample)
    metadata = pickle.load(open(samples.loc[selector, "pickle"].values[0], "rb"))
    atom_plddts = metadata["atom_plddts"]
    cif = samples.loc[selector, "cif"].values[0]
    cif_plddt = str(Path(cif).with_stem(Path(cif).stem + "_plddt"))
    parser = CIF2CIF(cif, cif_plddt)
    structure = parser.read()
    structure.set_annotation(
        "b_factor", atom_plddts
    )  # 不能使用 structure.b_factor = atom_plddt
    parser.write(structure)
