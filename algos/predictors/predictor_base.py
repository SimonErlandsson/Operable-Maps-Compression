from abc import ABC, abstractmethod


class PredictorFunction(ABC):

    @abstractmethod
    def predict(self, *predictors):
        pass
