
from algos.predictors.predictor_base import PredictorFunction

class SimplePredictor(PredictorFunction):

    def predict(self, prev_coords):
        return sum(prev_coords) / len(prev_coords)
