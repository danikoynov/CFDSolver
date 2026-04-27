import math

class Airfoil:
    def __init__(self, 
                 code: str, 
                 chord: float):
        if len(code) != 4:
            raise ValueError("NACA code must be 4 digits")
        
        self.chord = chord
        self.code = code
        self.m, self.p, self.t = self.decode(code)

    @staticmethod
    def decode(code):   
        m = int(code[0]) / 100.0
        p = int(code[1]) / 10.0
        t = int(code[2: 4]) / 100.0
        
        print("Check decode")
        print(f"m: {m} | p: {p} | t: {t}")
        return m, p, t
        
        
    def half_thickness(self, X: float):
        """
        Returns the half thickness at the point with 
        coordinate normalized coordinate X
        """
        y = 5 * self.t * self.chord * (
            0.2969 * math.sqrt(X) - 
            0.1260 * X -
            0.3516 * (X ** 2) +
            0.2843 * (X ** 3) -
            0.1015 * (X ** 4)
        )
        
        return y
    
    def mean_camber_line(self, X):
        """
        Returns the mean camber at the point with 
        coordinate normalized coordinate X
        """
        if X < self.p:
            y_c = self.chord * self.m / (self.p ** 2) * (
                2 * self.p * X - X ** 2
            )
        else:
            y_c = self.chord * self.m / ((1 - self.p) ** 2) * (
                (1 - 2 * self.p) + 
                2 * self.p * X -
                X ** 2
            )

        return y_c

