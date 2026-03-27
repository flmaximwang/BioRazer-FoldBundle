import shutil, json, pickle, warnings
from dataclasses import dataclass, field
from pathlib import Path
import pandas as pd
import numpy as np
from biorazer.structure.io.protein import PDB2CIF, CIF2PDB
from biorazer.design import SingleTest, Entry, Library


def _load_unique_npz_array(npz_path: str | Path):
    with np.load(npz_path) as npz_data:
        keys = npz_data.files
        if len(keys) != 1:
            raise ValueError(
                f"Expected exactly one array in {npz_path}, found {len(keys)}: {keys}"
            )
        return npz_data[keys[0]]


@dataclass
class Boltz2Sample(SingleTest):

    marker: str
    metadata_files: dict = field(
        default_factory=lambda: {
            "pickle": None,
            "cif": None,
            "pdb": None,
        }
    )  # Store any additional metadata related to the test
    metadata_file_suffixes: dict = field(
        default_factory=lambda: {
            "pickle": "_data.pkl",
            "cif": "_model.cif",
            "pdb": "_model.pdb",
        }
    )  # Expected suffixes for the metadata files, used for auto-discovery
    _pickle_cache: dict | None = None  # Cache for the loaded pickle
    simple_metrics: dict = field(
        default_factory=lambda: {
            "confidence_score": ["confidence_score"],
            "iptm": ["iptm"],
            "protein_iptim": ["protein_iptm"],
            "ligand_iptm": ["ligand_iptm"],
            "complex_plddt": ["complex_plddt"],
            "complex_iplddt": ["complex_iplddt"],
            "complex_pde": ["complex_pde"],
            "complex_ipde": ["complex_ipde"],
            "num_recycles": ["num_recycles"],
            "ptm": ["ptm"],
            "ranking_score": ["ranking_score"],
        }
    )


@dataclass
class Boltz2Result(Entry):

    # Properties inherited from Entry
    marker: str
    dir_path: Path | None = None
    tests: list[Boltz2Sample] | None = None
    test_type: type = Boltz2Sample
    dataframe: pd.DataFrame | None = None

    @classmethod
    def _guess_marker_from_dir(cls, dir_path):
        return super()._guess_marker_from_dir(dir_path)

    @classmethod
    def format(cls, src_dir: str | Path, target_dir: str | Path, keep_marker=True):
        src_dir, target_dir, old_marker, new_marker = cls._prepare_format_paths(
            src_dir, target_dir, keep_marker
        )

        models = list(src_dir.glob(f"{old_marker}_model_*"))
        n_samples = len(models)
        if n_samples == 0:
            raise ValueError(
                f"No model files found in {src_dir} with prefix {old_marker}_model_"
            )
        for i in range(n_samples):
            sample_dir = target_dir / f"{new_marker}_sample-{i}"
            if not sample_dir.exists():
                sample_dir.mkdir(parents=True, exist_ok=True)

            # Model files
            model_stem = f"{old_marker}_model_{i}"
            model_cif = src_dir / f"{model_stem}.cif"
            model_pdb = src_dir / f"{model_stem}.pdb"
            if model_cif.exists():
                model_raw_src = model_cif
                model_gen_src = model_pdb
                parser = CIF2PDB()
            elif model_pdb.exists():
                model_raw_src = model_pdb
                model_gen_src = model_cif
                parser = PDB2CIF()
            else:
                raise ValueError(
                    f"No model file found for sample {i} in {src_dir} with prefix {model_stem}"
                )
            model_raw_target = (
                sample_dir / f"{new_marker}_sample-{i}_model{model_raw_src.suffix}"
            )
            model_gen_target = (
                sample_dir / f"{new_marker}_sample-{i}_model{model_gen_src.suffix}"
            )
            shutil.copy(model_raw_src, model_raw_target)
            parser.input_file = model_raw_target
            parser.output_file = model_gen_target
            if isinstance(parser, CIF2PDB):
                structure = parser.read(extra_fields=["B_iso_or_equiv"])
                structure.set_annotation(
                    "b_factor", structure.B_iso_or_equiv.astype(float)
                )
            else:
                structure = parser.read()
            parser.write(structure)

            # Data files
            confidense_json = src_dir / f"confidence_{old_marker}_model_{i}.json"
            pae_npz = src_dir / f"pae_{old_marker}_model_{i}.npz"
            pde_npz = src_dir / f"pde_{old_marker}_model_{i}.npz"
            plddt_npz = src_dir / f"plddt_{old_marker}_model_{i}.npz"
            for file in [confidense_json, pae_npz, pde_npz, plddt_npz]:
                if not file.exists():
                    raise ValueError(
                        f"Expected data file {file} not found for sample {i} in {src_dir}"
                    )
            data_pickle_target = sample_dir / f"{new_marker}_sample-{i}_data.pkl"
            data = {}
            data.update(json.load(open(confidense_json, "r")))
            data["pae"] = _load_unique_npz_array(pae_npz)
            data["pde"] = _load_unique_npz_array(pde_npz)
            data["plddt"] = _load_unique_npz_array(plddt_npz)
            with open(data_pickle_target, "wb") as f:
                pickle.dump(data, f)
