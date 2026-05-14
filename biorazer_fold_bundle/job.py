from abc import abstractmethod
from dataclasses import dataclass, field


@dataclass
class Job:
    name: str

    @property
    @abstractmethod
    def dict(self):
        """
        Return a dictionary representation of the job, which can be used for serialization or execution.
        """


@dataclass
class JobBatch:

    jobs: list[Job]

    @abstractmethod
    def generate_requests(self):
        """
        Write the batch of jobs to disk or a database, making them ready for execution.
        """

    @abstractmethod
    def generate_command(self, *args, **kwargs):
        """
        Generate the command(s) needed to execute the batch of jobs, which can be run in a terminal or a job scheduler.
        """
