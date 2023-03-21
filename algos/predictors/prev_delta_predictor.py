from algos.predictors.predictor_base import PredictorFunction

class PrevDeltaPredictor(PredictorFunction):

    def predict(self, prev_coords):
        if len(prev_coords) == 0:
            return 0
        elif len(prev_coords) == 1:
            return prev_coords[0]
        else:
            return prev_coords[-1] - (prev_coords[-1] - prev_coords[-2])
