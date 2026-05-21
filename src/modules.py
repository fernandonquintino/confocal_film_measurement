import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
from scipy.optimize import curve_fit

def fit_func(x, a, b):
    return a*x**2 + b*x

def calib_pred(X, x0, x1, x2):
    X_pred = x0 + x1*X + x2*X**2
    return X_pred

class Calibration:
    """
    Class used for calibration.
    ...
    Attributes
    ----------
    X : pd.Series
        a formatted string to print out what the animal says
    y : pd.Series
        the name of the animal
        
    Methods
    -------
    models()
        Dictionary with each model as np.poly1d
    predictions()
        DataFrame with the predictions of each model
    errors(percentual=False)
        DataFrame with the relative error (NaN if the actual value is 0). If percentual is set to True, values are presented as percentages
    summary(percentual=False)
        DataFrame with both predictions and errors for each model
    """
    def __init__(self, X, y):
        self.X = X
        self.y = y
    
    def models(self):
        reg = LinearRegression().fit(self.X.values.reshape(-1, 1), self.y)
        reg2 = LinearRegression(fit_intercept=False).fit(self.X.values.reshape(-1, 1), self.y)
        reg3 = np.poly1d(np.polyfit(self.X, self.y, 2))
        reg4 = curve_fit(fit_func, self.X, self.y)

        Coef = {'1_i': [0, reg.coef_[0], reg.intercept_],
                '1_ni': [0, reg2.coef_[0], reg2.intercept_], 
                '2_i': [reg3[2], reg3[1], reg3[0]],
                '2_ni': [reg4[0][0], reg4[0][1], 0]}

        Models = {}
        for i in Coef.keys():
            Models[i] = np.poly1d(Coef[i])
            
        return Models
    
    def predictions(self):
        models = self.models()
        df = pd.DataFrame()
        for i in models.keys():
            df = pd.concat([df, pd.DataFrame({i: models[i](self.X)})], axis=1)
        df = df.add_prefix('pred_').round(decimals=3)
        return df
    
    def errors(self, percentual=False):
        y_corr = self.y
        df_pred = self.predictions()
        df_pred.columns = df_pred.columns.str.lstrip('pred_')
        error_df = df_pred.add_prefix('error_')
        for i in error_df.columns:
            error_df[i] = (self.y-error_df[i])/self.y
        error_df.replace([np.inf, -np.inf], np.nan, inplace=True)
        if percentual == True:
            error_df = error_df*100
            error_df = error_df.round(decimals=2)
        else:
            error_df = error_df.round(decimals=4)
        return error_df
    
    def summary(self, percentual=False):
        df_summary = pd.concat([self.predictions(), self.errors(percentual)], axis=1)
        return df_summary