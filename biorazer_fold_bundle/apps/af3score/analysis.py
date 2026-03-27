import re, json, shutil, pickle, warnings
from pathlib import Path
import pandas as pd
from tqdm import tqdm
from biorazer.structure.io import CIF2PDB
from ..analysis_basic import SinglePrediction, BatchPrediction
from ..alphafold3_ex.analysis import BatchAF3Local

warnings.filterwarnings(
    "ignore", category=UserWarning, module="biotite.structure.io.pdbx"
)


class SingleAF3Score(SinglePrediction):
    """
    This class is used to analyze the AF3 local prediction results.
    """

    def format_output(self, target_dir):

        target_marker_dir, _ = super().format_output(target_dir)

        for i, sample in enumerate(self.marker_dir_path.glob("seed-*_sample-*")):

            if not sample.is_dir():
                continue

            target_sample_dir_path = Path(target_marker_dir) / sample.stem
            if not target_sample_dir_path.exists():
                target_sample_dir_path.mkdir(parents=True)

            cif_src = sample / "model.cif"
            cif_target = target_sample_dir_path / (
                self.marker + "_" + sample.stem + "_model.cif"
            )
            shutil.copyfile(cif_src, cif_target)

            pdb_target = cif_target.with_suffix(".pdb")
            CIF2PDB(cif_src, pdb_target).convert()

            confidences_json_src = sample / "confidences.json"
            confidence_data = json.load(open(confidences_json_src, "r"))
            summary_confidences_json_src = sample / "summary_confidences.json"
            summary_confidences_data = json.load(
                open(summary_confidences_json_src, "r")
            )

            meta_data = {}
            for key in summary_confidences_data:
                meta_data[key] = summary_confidences_data[key]
            for key in confidence_data:
                meta_data[key] = confidence_data[key]
            meta_data_pickle_target = target_sample_dir_path / (
                self.marker + "_" + sample.stem + "_data.pickle"
            )
            with open(meta_data_pickle_target, "wb") as f:
                pickle.dump(meta_data, f)

        return target_marker_dir, SingleAF3Score(target_marker_dir)


class BatchAF3Score(BatchAF3Local):

    prediction_class = SingleAF3Score
