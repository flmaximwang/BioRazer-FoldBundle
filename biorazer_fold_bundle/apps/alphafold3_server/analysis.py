import warnings, json, shutil, pickle
from dataclasses import dataclass, field
from pathlib import Path
import pandas as pd
from biorazer.design.basic import *
from biorazer.structure.io.protein import CIF2PDB


@dataclass
class AF3ServerSample(SingleTest):
    """
    This class is used to analyze the results of a single AF3 server prediction run, which may contain multiple samples with the same seed.
    """

    marker: str
    metadata_files: dict = field(
        default_factory=lambda: {
            "pickle": None,
            "cif": None,
            "pdb": None,
            "msas_dir": None,
            "templates_dir": None,
        }
    )  # Store any additional metadata related to the test
    metadata_file_suffixes: dict = field(
        default_factory=lambda: {
            "pickle": "_data.pkl",
            "cif": "_model.cif",
            "pdb": "_model.pdb",
            "msas_dir": "_msas",
            "templates_dir": "_templates",
        }
    )  # Expected suffixes for the metadata files, used for auto-discovery
    _pickle_cache: dict | None = None  # Cache for the loaded pickle
    simple_metrics: dict = field(
        default_factory=lambda: {
            "fraction_disordered": ["fraction_disordered"],
            "iptm": ["iptm"],
            "num_recycles": ["num_recycles"],
            "ptm": ["ptm"],
            "ranking_score": ["ranking_score"],
        }
    )


@dataclass
class AF3ServerResult(Entry):

    # Properties inherited from Entry
    marker: str
    dir_path: Path | None = None
    tests: list[AF3ServerSample] | None = None
    test_type: type = AF3ServerSample
    dataframe: pd.DataFrame | None = None

    @classmethod
    def _guess_marker_from_dir(cls, dir_path):
        dir_path = Path(dir_path)
        dir_path_stem = dir_path.stem
        if dir_path_stem.startswith("fold_"):
            return dir_path_stem.replace("fold_", "")
        elif dir_path_stem.startswith("folds_"):
            raise ValueError(
                f"The directory name {dir_path_stem} seems to indicate that it contains multiple fold predictions. Please specify the marker name explicitly or use a directory that contains only one fold prediction."
            )
        else:
            return dir_path_stem

    @classmethod
    def format(cls, src_dir: str | Path, target_dir: str | Path, keep_marker=True):
        src_dir, target_dir, old_marker, new_marker = cls._prepare_format_paths(
            src_dir, target_dir, keep_marker
        )

        job_request_json = src_dir / f"fold_{old_marker}_job_request.json"
        job_request = json.load(open(job_request_json, "r"))[0]
        seed = job_request["modelSeeds"][0]
        for i in range(5):

            prefix = f"{new_marker}_seed-{seed}_sample-{i+1}"
            sample_seed_dir = target_dir / prefix
            if not sample_seed_dir.exists():
                sample_seed_dir.mkdir(parents=True)

            cif_src = src_dir / f"fold_{old_marker}_model_{i}.cif"
            cif_target = sample_seed_dir / f"{prefix}_model.cif"
            shutil.copyfile(cif_src, cif_target)

            # Add B-factor annotation (plddt) to the structure and write it to a PDB file
            pdb_target = cif_target.with_suffix(".pdb")
            parser = CIF2PDB(input_file=cif_src, output_file=pdb_target)
            structure = parser.read(extra_fields=["B_iso_or_equiv"])
            structure.set_annotation("b_factor", structure.B_iso_or_equiv.astype(float))
            parser.write(structure)

            pickle_data = {}
            pickle_data.update(job_request)
            full_data_json = src_dir / f"fold_{old_marker}_full_data_{i}.json"
            full_data = json.load(open(full_data_json, "r"))
            summary_confidences_json = (
                src_dir / f"fold_{old_marker}_summary_confidences_{i}.json"
            )
            summary_confidences_data = json.load(open(summary_confidences_json, "r"))
            pickle_data.update(full_data)
            pickle_data.update(summary_confidences_data)
            metadata_pickle_target = sample_seed_dir / f"{prefix}_data.pkl"
            with open(metadata_pickle_target, "wb") as f:
                pickle.dump(pickle_data, f)

            msas_dir_src = src_dir / f"msas"
            if not msas_dir_src.exists():
                warnings.warn(
                    f"MSAs directory {msas_dir_src} does not exist. Skipping copying MSAs directory for sample {prefix}."
                )
            else:
                msas_dir_target = sample_seed_dir / f"{prefix}_msas"
                if not msas_dir_target.exists():
                    msas_dir_target.mkdir(parents=True)
                for item in msas_dir_src.iterdir():
                    if item.is_file():
                        shutil.copyfile(item, msas_dir_target / item.name)
            templates_dir_src = src_dir / f"templates"
            if not templates_dir_src.exists():
                warnings.warn(
                    f"Templates directory {templates_dir_src} does not exist. Skipping copying templates directory for sample {prefix}."
                )
            else:
                templates_dir_target = sample_seed_dir / f"{prefix}_templates"
                if not templates_dir_target.exists():
                    templates_dir_target.mkdir(parents=True)
                for item in templates_dir_src.iterdir():
                    if item.is_file():
                        shutil.copyfile(item, templates_dir_target / item.name)
