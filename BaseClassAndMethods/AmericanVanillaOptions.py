from BaseClassAndMethods.BlackScholesBaseFormula import GeneralBlackScholesMertonFormula
from abc import ABC, abstractmethod
from math import log, sqrt, exp
from scipy.stats import norm

class BAWAmericanOptionsApprox(GeneralBlackScholesMertonFormula):
    """
    BaroneAdesiAndWhaley
    http://www.deriscope.com/docs/Barone_Adesi_Whaley_1987.pdf
    """
    def __init__(self, S: float, X: float, b:float, r: float, T: float, sigma: float, call_put: str) -> float:
        super().__init__(S, X, b, r, T, sigma, call_put)

    def Kc(self):
        """
        Newton-Raphson algorithm for Critical Comdty Price of a call option
        """
        

    @property
    def price(self):
        if self.b >= self.r:
            return super().price
        else:

        