import os, yaml, shutil
from dataclasses import dataclass, field
from pathlib import Path
from biorazer_fold_bundle.job import Job, JobBatch


class Boltz2Modification:

    def __init__(self, position: int, ccd: str):
        self.position = position
        self.ccd = ccd

    @property
    def dict(self):
        return {
            "position": self.position,
            "ccd": self.ccd,
        }


@dataclass
class Boltz2Sequence:
    """
    Parameters
    ----------
    entity_type: str
        The type of the sequence. Must be one of "protein", "dna", "rna", or "ligand".
    id: list of str
        The identifier(s) for the sequence.
    sequence: str, optional
        The amino acid or nucleotide sequence. Required for "protein", "dna", and "rna" entity types.
    msa: str, optional
        The multiple sequence alignment in a3m format. Optional for "protein", "dna", and "rna" entity types.
    smiles: str, optional
        The SMILES representation of the ligand. Required for "ligand" entity type if ccd is not provided.
    ccd: str, optional
        The CCD code of the ligand. Required for "ligand" entity type if smiles is not provided.
    modifications: list of Boltz2Modification, optional
        A list of modifications for the sequence. Applicable for "protein", "dna", and "rna" entity types.
    """

    entity_type: str
    id: list[str]
    sequence: str = None
    msa: str = None
    smiles: str = None
    ccd: str = None
    modifications: list[Boltz2Modification] = field(default_factory=list)

    def check_integrity(self):
        if self.entity_type == "protein":
            if self.sequence is None:
                raise ValueError("Sequence is required for protein entity type.")
            if self.msa is None:
                raise ValueError("MSA is required for protein entity type.")
        elif self.entity_type in ["dna", "rna"]:
            if self.sequence is None:
                raise ValueError("Sequence is required for dna and rna entity types.")
        elif self.entity_type == "ligand":
            if self.smiles is None and self.ccd is None:
                raise ValueError(
                    "Either SMILES or CCD code is required for ligand entity type."
                )
        else:
            raise ValueError(
                "Invalid entity type. Entity type must be 'protein', 'dna', 'rna', or 'ligand'."
            )

    def check_msa_type(self):
        if self.entity_type != "protein":
            raise ValueError("MSA is only applicable for protein entity type.")
        if self.msa is None:
            raise ValueError("MSA is not set.")
        if isinstance(self.msa, Path):
            self.msa = str(self.msa)
        if Path(self.msa).suffix != ".a3m":
            return "string"
        else:
            return "file"

    def seq_from_msa(self):
        msa_type = self.check_msa_type()
        if msa_type == "file":
            self.read_msa()
        flag = False
        for line in self.msa.splitlines():
            if (
                not line.startswith(">")
                and not line.startswith("#")
                and len(line.strip()) > 0
            ):
                self.sequence = line.strip()
                flag = True
                break
        if not flag:
            raise ValueError("No sequence found in MSA.")

    def msa_from_seq(self):
        """Actually no MSA, just write the sequence in a3m format."""
        self.msa = f">query\n{self.sequence}\n"

    def read_msa(self):
        msa_type = self.check_msa_type()
        if msa_type != "file":
            raise ValueError("MSA must be provided as a file for protein entity type.")
        with open(self.msa, "r") as f:
            self.msa = f.read()

    def write_msa(self, a3m: str | Path):
        msa_type = self.check_msa_type()
        if msa_type != "string":
            raise ValueError(
                "MSA must be provided as a string for protein entity type."
            )
        with open(a3m, "w") as f:
            f.write(self.msa)
        self.msa = str(a3m)

    @property
    def dict(self):
        self.check_integrity()
        output_dict = {
            "id": self.id,
        }
        if self.entity_type == "protein":
            if self.check_msa_type() != "file":
                raise ValueError(
                    "MSA must be provided as a file for protein entity type when exporting."
                )
            output_dict = {
                self.entity_type: {
                    "id": self.id,
                    "sequence": self.sequence,
                }
            }
            output_dict[self.entity_type]["msa"] = self.msa
            if len(self.modifications) > 0:
                output_dict[self.entity_type]["modifications"] = [
                    mod.dict() for mod in self.modifications
                ]
            return output_dict
        elif self.entity_type in ["dna", "rna"]:
            if self.sequence is None:
                raise ValueError("Sequence is required for dna and rna entity types.")
            output_dict = {
                self.entity_type: {
                    "id": self.id,
                    "sequence": self.sequence,
                }
            }
            if len(self.modifications) > 0:
                output_dict[self.entity_type]["modifications"] = [
                    mod.dict() for mod in self.modifications
                ]
            return output_dict
        elif self.entity_type == "ligand":
            # 优先使用 SMILES
            if not self.smiles is None:
                return {
                    self.entity_type: {
                        "id": self.id,
                        "smiles": self.smiles,
                    }
                }
            else:
                return {
                    self.entity_type: {
                        "id": self.id,
                        "ccd": self.ccd,
                    }
                }
        else:
            raise ValueError(
                "Invalid entity type. Entity type must be 'protein', 'dna', 'rna', or 'ligand'"
            )


@dataclass
class Boltz2Constraint:
    pass


@dataclass
class Boltz2BondAtom:
    chain_id: str
    res_id: int
    atom_name: str

    @property
    def dict(self):
        return [self.chain_id, self.res_id, self.atom_name]


@dataclass
class Boltz2ConstraintBond(Boltz2Constraint):
    atom1: Boltz2BondAtom
    atom2: Boltz2BondAtom

    @property
    def dict(self):
        return {"bond": {"atom1": self.atom1.dict, "atom2": self.atom2.dict}}


@dataclass
class Boltz2Contact:
    chain_id: str
    res_id: int | None = None
    atom_name: str | None = None

    def is_legal(self):
        if self.res_id is None and self.atom_name is None:
            raise ValueError("At least one of res_id and atom_name must be provided.")
        if self.res_id is not None and self.atom_name is not None:
            raise ValueError("Only one of res_id and atom_name can be provided.")

    @property
    def dict(self):
        self.is_legal()
        if self.res_id is not None:
            return [self.chain_id, self.res_id]
        else:
            return [self.chain_id, self.atom_name]


@dataclass
class Boltz2ConstraintPocket(Boltz2Constraint):
    binder: str  # Chain id
    contacts: list[Boltz2Contact]
    max_distance: float | None = None
    # max_distance specifies the maximum distance (in Angstrom, supported between
    # 4A and 20A with 6A as default) between any atom in the binder and any atom
    # in each of the contacts elements.
    force: bool = False
    # if force is set to true (default is false), a potential
    # will be used to enforce the contact constraint

    @property
    def dict(self):
        return {
            "pocket": {
                "binder": self.binder,
                "contacts": [contact.dict for contact in self.contacts],
                "max_distance": self.max_distance,
                "force": self.force,
            }
        }


@dataclass
class Boltz2ConstraintContact(Boltz2Constraint):
    token1: Boltz2Contact
    token2: Boltz2Contact
    max_distance: float | None = None
    # max_distance specifies the maximum distance (in Angstrom, supported between
    # 4A and 20A with 6A as default) between any atom in the binder and any atom
    # in each of the contacts elements.
    force: bool = False
    # if force is set to true (default is false), a potential
    # will be used to enforce the contact constraint

    @property
    def dict(self):
        return {
            "contact": {
                "token1": self.token1.dict,
                "token2": self.token2.dict,
                "max_distance": self.max_distance,
                "force": self.force,
            }
        }


@dataclass
class Boltz2Template:
    cif: str | None = None  # The path to the CIF file of the template structure.
    pdb: str | None = None  # The path to the PDB file of the template structure.
    # Specify which chain(s) need template information.
    chain_id: str | list[str] | None = None
    # Specify which chain(s) in the template structure correspond to the target sequence(s).
    template_id: str | list[str] | None = None
    force: bool | None = None
    threshold: float | None = None

    def __post_init__(self):
        if self.cif is None and self.pdb is None:
            raise ValueError("Either cif or pdb must be provided for template.")
        if self.cif is not None and self.pdb is not None:
            raise ValueError("Only one of cif or pdb can be provided for template.")
        if self.template_id is not None and self.chain_id is None:
            raise ValueError(
                "If template_id is provided, chain_id must also be provided."
            )

    def dict(self):
        output_dict = {}
        if self.cif is not None:
            output_dict["cif"] = self.cif
        if self.pdb is not None:
            output_dict["pdb"] = self.pdb
        if self.chain_id is not None:
            output_dict["chain_id"] = self.chain_id
        if self.template_id is not None:
            output_dict["template_id"] = self.template_id
        if self.force is not None:
            output_dict["force"] = self.force
        if self.threshold is not None:
            output_dict["threshold"] = self.threshold
        return output_dict


@dataclass
class Boltz2Job(Job):

    name: str
    sequences: list[Boltz2Sequence] | None = None
    constraints: list[Boltz2Constraint] | None = None
    templates: list[Boltz2Template] | None = None

    @property
    def dict(self):
        output_dict = {
            "sequences": [],
        }
        for seq in self.sequences:
            output_dict["sequences"].append(seq.dict)
        if self.constraints and len(self.constraints) > 0:
            output_dict["constraints"] = []
            for con in self.constraints:
                output_dict["constraints"].append(con.dict)
        if self.templates and len(self.templates) > 0:
            output_dict["templates"] = []
            for template in self.templates:
                output_dict["templates"].append(template.dict)
        return output_dict


@dataclass
class Boltz2JobBatch(JobBatch):

    jobs: list[Boltz2Job] | None = None

    def generate_requests(self, target_dir):
        target_dir_path = Path(target_dir)
        if not target_dir_path.exists():
            target_dir_path.mkdir(parents=True)

        ori_cwd = os.getcwd()
        os.chdir(target_dir)
        msa_dir = Path("msas")
        request_dir = Path("requests")
        for c_dir in [msa_dir, request_dir]:
            if not c_dir.exists():
                c_dir.mkdir(parents=True)
        for job in self.jobs:
            job_name = job.name
            for sequence in job.sequences:
                if sequence.entity_type == "protein":
                    seq_msa_dir = msa_dir / f"{job_name}"
                    if not seq_msa_dir.exists():
                        seq_msa_dir.mkdir(parents=True)
                    msa_type = sequence.check_msa_type()
                    if msa_type == "string":
                        a3m = seq_msa_dir / f"{job.name}_{''.join(sequence.id)}.a3m"
                        sequence.write_msa(a3m)
                    elif msa_type == "file":
                        if Path(sequence.msa).exists():
                            target_a3m = (
                                seq_msa_dir / f"{job.name}_{''.join(sequence.id)}.a3m"
                            )
                            shutil.copy(sequence.msa, target_a3m)
                            sequence.msa = str(target_a3m)
                        else:
                            raise FileNotFoundError(
                                f"MSA file {sequence.msa} not found."
                            )
                    else:
                        raise ValueError(
                            f"Invalid MSA type {msa_type}. MSA type must be 'string' or 'file'."
                        )
            output_dict: dict = job.dict
            with open(request_dir / f"{job.name}.yaml", "w") as f:
                yaml.dump(output_dict, f, indent=2)
        os.chdir(ori_cwd)

    def generate_command(
        self,
        input_file: str | Path,
        output_dir: str | Path = ".",
        recycling_steps: int | None = None,
        sampling_steps: int | None = None,
        diffusion_samples: int | None = 5,
        output_format: str | None = None,
        cache_dir: str | Path = "~/.boltz",
        use_msa_server=False,
        msa_server_url: str = "https://api.colabfold.com",
    ):
        """
        See https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md for details

        Parameters
        ----------
        input_file: str or Path
            The path to the input YAML/fasta file containing the job request.
        output_dir: str or Path, optional
            The directory where the output files will be saved. Default is the current directory.
        cache_dir: str or Path, optional
            The directory where the cache files will be stored. Default is "~/.boltz".
        use_msa_server: bool, optional
            Whether to use the MSA server. Default is False.
        msa_server_url: str, optional
            The URL of the MSA server. Default is "https://api.colabfold.com".
        output_format: str, optional
            The format of the output files. Default is "mmcif".
        recycling_steps: int, optional
            The number of recycling steps. Default is 3.
        sampling_steps: int, optional
            The number of sampling steps. Default is 200.


        """
        commands = []
        commands.extend(["boltz", "predict"])
        commands.append(f'"{input_file}"')
        commands.append(f'--out_dir "{output_dir}"')
        commands.append(f'--cache "{cache_dir}"')
        if use_msa_server:
            commands.append("--use-msa-server")
            commands.append(f'--msa-server-url "{msa_server_url}"')
        if output_format:
            commands.append(f'--output_format "{output_format}"')
        if recycling_steps is not None:
            commands.append(f"--recycling_steps {recycling_steps}")
        if sampling_steps is not None:
            commands.append(f"--sampling_steps {sampling_steps}")
        if diffusion_samples is not None:
            commands.append(f"--diffusion_samples {diffusion_samples}")

        return " ".join(commands)
