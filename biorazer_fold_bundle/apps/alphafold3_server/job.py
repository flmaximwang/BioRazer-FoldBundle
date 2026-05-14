import json
from abc import abstractmethod
from dataclasses import dataclass, field
from pathlib import Path

from biorazer_fold_bundle.job import Job, JobBatch


@dataclass
class AF3ServerSequence:
    count: int = 1

    @property
    @abstractmethod
    def dict(self):
        """
        Return a dictionary representation of a sequence entry accepted by the
        AlphaFold3 server request format.
        """


@dataclass
class AF3ServerProteinChain(AF3ServerSequence):
    sequence: str = ""
    count: int = 1

    @property
    def dict(self):
        if not self.sequence:
            raise ValueError(
                "Sequence is required for an AlphaFold3 server protein chain."
            )
        return {
            "proteinChain": {
                "sequence": self.sequence,
                "count": self.count,
            }
        }


@dataclass
class AF3ServerJob(Job):
    name: str
    sequences: list[AF3ServerSequence] = field(default_factory=list)
    model_seeds: list[int] = field(default_factory=list)

    @property
    def dict(self):
        if not self.name:
            raise ValueError("Job name is required.")
        if len(self.sequences) == 0:
            raise ValueError(
                "At least one sequence is required for an AlphaFold3 server job."
            )
        return {
            "name": self.name,
            "modelSeeds": self.model_seeds,
            "sequences": [sequence.dict for sequence in self.sequences],
        }


@dataclass
class AF3ServerJobBatch(JobBatch):
    jobs: list[AF3ServerJob] | None = None

    def generate_requests(
        self,
        target_dir: str | Path,
        request_num_per_json: int = 30,
        target_prefix: str = "af3_server_requests",
    ):
        if request_num_per_json < 1:
            raise ValueError("request_num_per_json must be at least 1.")
        if not target_prefix:
            raise ValueError("target_prefix must not be empty.")

        target_dir = Path(target_dir)
        target_dir.mkdir(parents=True, exist_ok=True)

        request_dir = target_dir / "requests"
        request_dir.mkdir(parents=True, exist_ok=True)

        jobs = self.jobs or []
        for start_idx in range(0, len(jobs), request_num_per_json):
            chunk = jobs[start_idx : start_idx + request_num_per_json]
            payload = [job.dict for job in chunk]
            file_idx = start_idx // request_num_per_json + 1
            file_name = f"{target_prefix}_{file_idx}.json"

            with open(request_dir / file_name, "w") as handle:
                json.dump(payload, handle, indent=4)

    def generate_command(self, *args, **kwargs):
        raise NotImplementedError(
            "AlphaFold3 Server requests are intended for manual upload; no CLI command is available."
        )
