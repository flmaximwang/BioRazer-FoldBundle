import json, os, shutil
from pathlib import Path
from abc import abstractmethod
from dataclasses import dataclass, field
from .._basic import SingleExecution, BatchExecution


@dataclass
class AF3Glycan:
    residues: str
    position: int

    @property
    def json(self):
        return {"residues": self.residues, "position": self.position}


@dataclass
class AF3Modification:
    ptm_type: str
    ptm_position: int

    @property
    def json(self):
        return {"ptmType": self.ptm_type, "ptmPosition": self.ptm_position}


@dataclass
class AF3Template:
    mmcif: str
    query_indices: list[int]
    template_indices: list[int]

    @staticmethod
    def from_mmcif(mmcif: str, query_indices: list[int], template_indices: list[int]):
        with open(mmcif, "r") as f:
            return AF3Template(f.read(), query_indices, template_indices)

    @property
    def json(self):
        return {
            "mmcif": self.mmcif,
            "queryIndices": self.query_indices,
            "templateIndices": self.template_indices,
        }


@dataclass
class AF3Sequence:
    id: list[str] | str

    @property
    @abstractmethod
    def json(self):
        """
        Implement in subclasses to return a dict that can easily be parsed by json, yaml, etc.
        """


@dataclass
class AF3ProteinLocal(AF3Sequence):
    """
    Parameters
    ----------
    id : list[str] | str
        Identifier(s) for the protein sequence.
    sequence : str
        The amino acid sequence of the protein.
    unpaired_msa : str, optional
        The unpaired multiple sequence alignment (MSA) in string format. Defaults to an empty string.
    unpaired_msa_path : str, optional
        The file path to the unpaired MSA. Defaults to an empty string.
    paired_msa : str, optional
        The paired multiple sequence alignment (MSA) in string format. Defaults to an empty string.
    paired_msa_path : str, optional
        The file path to the paired MSA. Defaults to an empty string.
    modifications : list[AF3Modification], optional
        A list of post-translational modifications for the protein. Defaults to an empty list.
    templates : list[AF3Template], optional
        A list of structural templates for the protein. Defaults to an empty list.
    """

    id: list[str] | str
    sequence: str
    unpaired_msa: str = ""
    unpaired_msa_path: str = ""
    paired_msa: str = ""
    paired_msa_path: str = ""
    modifications: list[AF3Modification] = field(default_factory=lambda: [])
    templates: list[AF3Template] = field(default_factory=lambda: [])

    def set_no_paired_msa(self):
        self.paired_msa = f">query\n{self.sequence}"

    def set_no_unpaired_msa(self):
        self.unpaired_msa = f">query\n{self.sequence}"

    def set_no_msa(self):
        self.set_no_paired_msa()
        self.set_no_unpaired_msa()

    def set_unpaired_msa(self, a3m: str):
        with open(a3m, "r") as f:
            self.unpaired_msa = f.read()

    def set_paired_msa(self, a3m: str):
        with open(a3m, "r") as f:
            self.paired_msa = f.read()

    def seq_from_unpaired_msa(self, a3m: str):
        with open(a3m, "r") as f:
            for line in f:
                if not line.startswith(">"):
                    self.sequence = line.strip()
                    break
        self.unpaired_msa_path = a3m

    def seq_from_paired_msa(self, a3m: str):
        with open(a3m, "r") as f:
            for line in f:
                if not line.startswith(">"):
                    self.sequence = line.strip()
                    break
        self.paired_msa_path = a3m

    def export_unpaired_msa(self, a3m: str):
        with open(a3m, "w") as f:
            f.write(self.unpaired_msa)
        self.unpaired_msa_path = a3m
        self.unpaired_msa = ""

    def export_paired_msa(self, a3m: str):
        with open(a3m, "w") as f:
            f.write(self.paired_msa)
        self.paired_msa_path = a3m
        self.paired_msa = ""

    def export_msa(self, unpaired_a3m: str, paired_a3m: str):
        self.export_unpaired_msa(unpaired_a3m)
        self.export_paired_msa(paired_a3m)

    @property
    def json(self):
        if self.paired_msa:
            raise ValueError(
                "Paired MSA string is set. Please export it to a file and clear the string before output."
            )
        if self.unpaired_msa:
            raise ValueError(
                "Unpaired MSA string is set. Please export it to a file and clear the string before output."
            )
        return {
            "protein": {
                "id": self.id,
                "sequence": self.sequence,
                "unpairedMsaPath": self.unpaired_msa_path,
                "pairedMsaPath": self.paired_msa_path,
                "modifications": [mod.json for mod in self.modifications],
                "templates": [temp.json for temp in self.templates],
            },
        }


@dataclass
class AF3Atom:
    chain_id: str
    res_id: int
    atom_name: str

    @property
    def json(self):
        return [self.chain_id, self.res_id, self.atom_name]


@dataclass
class AF3BondedAtomPair:
    atom1: AF3Atom
    atom2: AF3Atom

    @property
    def json(self):
        return [self.atom1.json, self.atom2.json]


@dataclass
class SingleAF3LocalEx(SingleExecution):

    name: str
    modelSeeds: list[int] = field(default_factory=lambda: [42])
    sequences: list[AF3Sequence] = field(default_factory=lambda: [])
    bonded_atom_pairs: list[AF3BondedAtomPair] = field(default_factory=lambda: [])
    userCCD: str = ""
    dialect: str = "alphafold3"  # Required
    version: int = 3  # Required

    @property
    def json(self):
        res = {
            "name": self.name,
            "modelSeeds": self.modelSeeds,
            "sequences": [seq.json for seq in self.sequences],
            "bondedAtomPairs": [pair.json for pair in self.bonded_atom_pairs],
            "dialect": self.dialect,
            "version": self.version,
        }
        if self.userCCD:
            res["userCCD"] = self.userCCD
        return res


@dataclass
class BatchAF3LocalEx(BatchExecution):
    executions: list[SingleAF3LocalEx]

    def write(self, target_dir):
        target_dir_path = Path(target_dir)
        if not target_dir_path.exists():
            target_dir_path.mkdir(parents=True, exist_ok=True)
        ori_cwd = os.getcwd()
        os.chdir(target_dir)

        for exe in self.executions:
            for seq in exe.sequences:
                if isinstance(seq, AF3ProteinLocal):
                    if seq.unpaired_msa:
                        seq.export_unpaired_msa(f"./{exe.name}_{seq.id}_unpaired.a3m")
                    else:
                        if not Path(seq.unpaired_msa_path).is_relative_to(
                            target_dir_path
                        ):
                            shutil.copy(
                                seq.unpaired_msa_path,
                                target_dir_path / f"{exe.name}_{seq.id}_unpaired.a3m",
                            )
                            seq.unpaired_msa_path = (
                                f"./{exe.name}_{seq.id}_unpaired.a3m"
                            )
                    if seq.paired_msa:
                        seq.export_paired_msa(f"./{exe.name}_{seq.id}_paired.a3m")
                    else:
                        if not Path(seq.paired_msa_path).is_relative_to(
                            target_dir_path
                        ):
                            shutil.copy(
                                seq.paired_msa_path,
                                target_dir_path / f"{exe.name}_{seq.id}_paired.a3m",
                            )
                            seq.paired_msa_path = f"./{exe.name}_{seq.id}_paired.a3m"
            res = exe.json
            with open(target_dir_path / f"{exe.name}.json", "w") as f:
                json.dump(res, f, indent=4)

        os.chdir(ori_cwd)
