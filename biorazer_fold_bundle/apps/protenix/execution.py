import json, os, shutil, re
from pathlib import Path
from contextlib import contextmanager
from dataclasses import dataclass, field
from ..execution_basic import SingleExecution, BatchExecution


@contextmanager
def change_dir(target_dir: str | Path):
    ori_dir = os.getcwd()
    os.chdir(target_dir)
    try:
        yield
    finally:
        os.chdir(ori_dir)


@dataclass
class ProtenixSequence:
    pass


@dataclass
class ProtenixModification:
    ptm_type: str
    ptm_position: int

    def output(self):
        return dict(ptmType=self.ptm_type, ptmPosition=self.ptm_position)


@dataclass
class ProtenixProteinChain(ProtenixSequence):
    sequence: str
    count: int
    pairing_msa: str = ""
    non_pairing_msa: str = ""
    precomputed_msa_dir: str = field(default_factory=lambda: "")
    modifications: list[ProtenixModification] = field(default_factory=list)
    pairing_db: str = "uniref100"  # You shuold not change it

    def output(self):
        return {
            "proteinChain": {
                "sequence": self.sequence,
                "count": self.count,
                "msa": {
                    "precomputed_msa_dir": self.precomputed_msa_dir,
                    "pairing_db": "uniref100",
                },
                "modifications": [mod.output() for mod in self.modifications],
            }
        }

    def set_pairing_msa(self, a3m: str):
        self.pairing_msa = a3m

    def set_non_pairing_msa(self, a3m: str):
        self.non_pairing_msa = a3m

    def set_msa_dir(self, msa_dir: str):
        self.precomputed_msa_dir = msa_dir

    def seq_from_pairing_msa(self, a3m: str):
        with open(a3m, "r") as f:
            for line in f:
                if not line.startswith(">"):
                    self.sequence = line.strip()
                    break
        self.set_pairing_msa(a3m)

    def seq_from_non_pairing_msa(self, a3m: str):
        with open(a3m, "r") as f:
            for line in f:
                if not line.startswith(">"):
                    self.sequence = line.strip()
                    break
        self.set_non_pairing_msa(a3m)

    def set_no_pairing_msa(self):
        self.pairing_msa = f">query\n{self.sequence}"

    def set_no_non_pairing_msa(self):
        self.non_pairing_msa = f">query\n{self.sequence}"

    def set_no_msa(self):
        self.set_no_pairing_msa()
        self.set_no_non_pairing_msa()

    def check_msa_content_type(self, msa_content):
        if re.match(r">[^\n]*query.*", msa_content):  # msa is file content
            return "content"
        elif Path(msa_content).exists():  # msa is valid file path
            return "file"
        else:
            raise ValueError(f"Invalid msa content: {msa_content[:200]}")

    def make_msa_absolute(self):
        if self.check_msa_content_type(self.pairing_msa) == "file":
            self.pairing_msa = str(Path(self.pairing_msa).absolute())
        if self.check_msa_content_type(self.non_pairing_msa) == "file":
            self.non_pairing_msa = str(Path(self.non_pairing_msa).absolute())

    def export_msa_content(self, msa_content, msa_type: str):
        """

        Parameters
        ----------
        msa_type: str
            "pairing" or "non_pairing"
        """
        msa_dir_path = Path(self.precomputed_msa_dir)
        if not msa_dir_path.exists():
            msa_dir_path.mkdir(parents=True)
        if re.match(r">[^\n]*query.*", msa_content):  # msa is file content
            with open(msa_dir_path / f"{msa_type}.a3m", "w") as f:
                f.write(msa_content)
        elif Path(msa_content).exists():  # msa is valid file path
            shutil.copy(msa_content, msa_dir_path / f"{msa_type}.a3m")
        else:
            raise ValueError(f"Invalid msa content: {msa_content[:200]}")

    def export_msa(self):
        if self.precomputed_msa_dir == "":
            raise ValueError("Please set precomputed_msa_dir before exporting msa")
        self.export_msa_content(self.pairing_msa, "pairing")
        self.pairing_msa = ""
        self.export_msa_content(self.non_pairing_msa, "non_pairing")
        self.non_pairing_msa = ""


@dataclass
class ProtenixColvalentBond:
    entity1: str
    copy1: int
    position1: str
    atom1: str
    entity2: str
    copy2: int
    position2: str
    atom2: str

    def output(self):
        return {
            "entity1": self.entity1,
            "copy1": self.copy1,
            "position1": self.position1,
            "atom1": self.atom1,
            "entity2": self.entity2,
            "copy2": self.copy2,
            "position2": self.position2,
            "atom2": self.atom2,
        }


@dataclass
class SingleProtenixEx(SingleExecution):
    name: str
    sequences: list[ProtenixSequence] = field(default_factory=list)
    covalent_bonds: list[ProtenixColvalentBond] = field(default_factory=list)

    def export_msa(self, target_dir: str | Path):
        """
        Export MSA for each protein chain to target_dir.
        """
        target_dir_path = Path(target_dir)
        if not target_dir_path.exists():
            target_dir_path.mkdir(parents=True)

        counter = -1
        for seq in self.sequences:
            if isinstance(seq, ProtenixProteinChain):
                counter += 1
                seq.make_msa_absolute()
                with change_dir(target_dir):
                    seq.set_msa_dir(f"./{self.name}_msa_{counter}")
                    seq.export_msa()

    def json(self):
        return dict(
            name=self.name,
            sequences=[seq.output() for seq in self.sequences],
            covalent_bonds=[bond.output() for bond in self.covalent_bonds],
        )


@dataclass
class BatchProtenixEx(BatchExecution):
    executions: list[SingleProtenixEx]

    def write(self, target_dir: str, num_per_file: int | None = None):
        target_dir_path = Path(target_dir)
        if not target_dir_path.exists():
            target_dir_path.mkdir(parents=True)
        for exe in self.executions:
            exe.export_msa(target_dir)
        res = [exe.json() for exe in self.executions]
        if num_per_file is None:
            num_per_file = len(res)
        with change_dir(target_dir):
            for i in range(0, len(res), num_per_file):
                with open(f"BatchProtenixEx_{i//num_per_file+1}.json", "w") as f:
                    json.dump(res[i : i + num_per_file], f, indent=4)
