import pickle
from abc import ABC, abstractmethod
from datetime import datetime
from pathlib import Path
from sys import getsizeof

import numpy as np

from pwmv_ecog.paths import ECoGPaths
from pwmv_ecog.plotting import Subjects


class Analysis(ABC):
    def __init__(self, subjects: Subjects, outputPath: Path) -> None:
        self.subjects = subjects
        self.rootFolder = outputPath

    def __sizeof__(self) -> int:
        # Return size of object in bytes
        sizeSum = super().__sizeof__() + getsizeof(self.subjects) + getsizeof(self.rootFolder)
        if hasattr(self, "data"):
            sizeSum += getsizeof(self.data)
        if hasattr(self, "results"):
            sizeSum += getsizeof(self.results)
        return sizeSum

    @property
    @abstractmethod
    def dataFolder(self) -> Path:
        pass

    @property
    @abstractmethod
    def saveFolder(self) -> Path:
        pass

    @property
    def name(self) -> str:
        return self.__class__.__name__

    @abstractmethod
    def run(self):
        # Execute analysis with current parameters
        pass

    @property
    def _saveFilename(self) -> str:
        return f"{self.name}_{self.make_detail_str()}.pkl"

    def save(self, path: Path = None):
        if not path:
            path = self.saveFolder / self._saveFilename
        with open(path, "wb") as f:
            pickle.dump(self, f)
        print(f"Saved {self.name} to {path}.")

    def load(self, pklPath: Path = None, runIfNotSaved=True):
        if not pklPath:
            pklPath = self.saveFolder / self._saveFilename
        if pklPath.is_file():
            rootFolder = self.rootFolder
            with open(pklPath, "rb") as f:
                self = pickle.load(f)
                self.rootFolder = rootFolder
            print(f"Loaded {self.name} from {pklPath}.")
            return self
        elif runIfNotSaved:
            print(f"No saved copy found at {pklPath}. Running analysis.")
            self.run()
            self.save(pklPath)
            return self
        else:
            raise FileNotFoundError(f"{pklPath} not found, and runIfNotSaved not True.")

    @abstractmethod
    def load_data(self):
        # Set .data and possibly other attributes
        pass

    @abstractmethod
    def make_detail_str(self) -> str:
        """
        Create a string that can be used to uniquely identify Analysis parameters.
        """
        # detail = "-".join([f"{val}" for val in kwargs.values()])
        # return detail
        pass

    @staticmethod
    def _basic_detail_str(**kwargs) -> str:
        """
        For simple make_detail_str implementations: joins all kwargs with a hyphen.
        """
        return "_".join([f"{val}" for val in kwargs.values()])

    @staticmethod
    @abstractmethod
    def _source_last_updated(sourceID) -> datetime | None:
        """
        Return a datetime object of the last time a data source was processed. None if not
        processed.

        Args:
            sourceID (Any): Identifier for data source to check.

        Returns:
            datetime | None: Datetime of last processing if available.
        """
        pass


class PWMVAnalysis(Analysis):
    # All data expected to be at 512 Hz
    dataSampleRate = 512

    def __init__(
        self,
        subjects: Subjects,
        outputPath: Path,
        smoothWindowSecs: float | None = 0.1,
    ) -> None:
        super().__init__(subjects, outputPath)
        self.smoothWindowSecs = smoothWindowSecs

    def __sizeof__(self) -> int:
        return (
            super().__sizeof__()
            + getsizeof(self.smoothWindowSecs)
            + getsizeof(self.dataSampleRate)
        )

    @staticmethod
    def _source_last_updated(patientID: int | str) -> datetime | None:
        """
        Return a datetime object of the last time the patient was processed. None if not processed.

        Args:
            patientID (int | str): Patient ID number.

        Returns:
            datetime | None: Datetime of last processing if available.
        """
        procPath = Path(ECoGPaths()._patient_path("path_processed_dir", f"Patient_{patientID}"))
        if (ymlPath := procPath / "processed_subinfo.yaml").is_file():
            return datetime.fromtimestamp(ymlPath.stat().st_mtime)
        else:
            return None

    @property
    def _update_before_date(self) -> datetime:
        if not hasattr(self, "_update_before"):
            sourceUpdates = [self._source_last_updated(subj) for subj in self.subjects]
            cleanDates = [dt for dt in sourceUpdates if dt is not None]
            if not cleanDates:
                self._update_before = datetime.max
            else:
                self._update_before = max(cleanDates)
        return self._update_before

    @property
    def _smoothingWindow(self) -> int | None:
        # Witdth of smoothing window in datapoints
        return (
            None
            if self.smoothWindowSecs is None
            else int(np.floor(self.smoothWindowSecs * self.dataSampleRate))
        )

    @property
    def dataFolder(self) -> Path:
        # Results data folder
        if not hasattr(self, "_dataFolder"):
            self._dataFolder = Path(self.rootFolder, "data")
            self._dataFolder.mkdir(parents=True, exist_ok=True)
        return self._dataFolder

    def load(self, pklPath: Path = None, runIfNotSaved=True):
        self = super().load(pklPath, runIfNotSaved)
        # Reset saved paths
        if hasattr(self, "_dataFolder"):
            delattr(self, "_dataFolder")
        return self
