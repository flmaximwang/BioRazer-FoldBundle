import re, json, shutil, pickle, warnings
from pathlib import Path
from dataclasses import dataclass, field
import pandas as pd
from tqdm import tqdm
from biorazer.structure.io.protein import CIF2PDB, CIF2CIF
from biorazer.design.basic.entry import *
from ...analysis import *


class SingleAF3Local(SinglePrediction):
    """
    This class is used to analyze the AF3 local prediction results.
    """

    def format_output(self, target_dir):

        target_marker_dir, _ = super().format_output(target_dir)
        entry_dir = self.dir_path

        for i, sample in enumerate(entry_dir.glob("seed-*_sample-*")):

            if not sample.is_dir():
                continue

            target_sample_dir_path = Path(target_marker_dir) / sample.stem
            if not target_sample_dir_path.exists():
                target_sample_dir_path.mkdir(parents=True)

            cif_src = sample / (entry_dir.stem + "_" + sample.stem + "_model.cif")
            cif_target = target_sample_dir_path / (
                entry_dir.stem + "_" + sample.stem + "_model.cif"
            )
            shutil.copyfile(cif_src, cif_target)

            pdb_target = cif_target.with_suffix(".pdb")
            parser = CIF2PDB(cif_src, pdb_target)
            structure = parser.read(extra_fields=["B_iso_or_equiv"])
            structure.set_annotation("b_factor", structure.B_iso_or_equiv.astype(float))
            parser.write(structure)

            confidences_json = sample / (
                entry_dir.stem + "_" + sample.stem + "_confidences.json"
            )
            confidence_data = json.load(open(confidences_json, "r"))
            summary_confidences_json = sample / (
                entry_dir.stem + "_" + sample.stem + "_summary_confidences.json"
            )
            summary_confidences_data = json.load(open(summary_confidences_json, "r"))

            meta_data = {}
            meta_data.update(summary_confidences_data)
            meta_data.update(confidence_data)
            meta_data_pickle_target = target_sample_dir_path / (
                entry_dir.stem + "_" + sample.stem + "_data.pickle"
            )
            with open(meta_data_pickle_target, "wb") as f:
                pickle.dump(meta_data, f)

        return target_marker_dir, SingleAF3Local(target_marker_dir)


class BatchAF3Local(BatchPrediction):

    single_pred_cls = SingleAF3Local
