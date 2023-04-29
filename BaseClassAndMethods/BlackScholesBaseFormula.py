from abc import ABC, abstractmethod

class BlackScholeNonDividend(ABC):
    """The most basic Black-Scholes Formula that values a European Option on a stock that does not pay dividends before the option's expiry date

    Args:
        S: Spot Price of the Stock
        X: Strike Price of the Option, usually denoted in textbooks as K
        r: Annualized Risk Free Interest Rate, usually the 10-Y Treasury Yield
        T: Time to Maturity in Years
        sigma: Annualized Historical Volatility (Standard Deviation) of the Stock Price 
        call_put: 'c' for call, 'p' for put
    
    Returns:
        Option Price for a non-dividend paying stock
        
    """

    def __init__(self, S: float, X: float, r: float, T: float, sigma: float, call_put: str) -> None:
        super().__init__()
        assert 