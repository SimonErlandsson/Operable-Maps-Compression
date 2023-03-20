from abc import ABC, abstractmethod


class PredictorFunction(ABC):

    #Function encapsulated to chunks
    @abstractmethod
    def predict(self, *predictors):
        pass
