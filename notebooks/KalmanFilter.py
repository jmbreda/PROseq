import numpy as np

class KalmanFilter(object):
    def __init__(self, F = None, B = None, H = None, Q = None, R = None, μ_0 = None, Σ_0 = None):

        if(F is None or H is None):
            raise ValueError("Set proper system dynamics.")
        
        # Parameters of the model
        self.m, self.n = H.shape # measurement and hidden state dimensions
        self.F = F # transition matrices
        self.H = H # Observation matrix
        self.B = np.zeros(1) if B is None else B # control matrix (optional)
        self.Q = np.eye(self.n) if Q is None else Q # transition noise covariance
        self.R = np.eye(self.m) if R is None else R # observation noise covariance
        self.S = None

        # save initial state for backward pass
        self.μ_0 = μ_0
        self.Σ_0 = Σ_0
        
        # predicted state
        self.μ_pred = None
        self.Σ_pred = None        

        # updated / current state
        self.μ = np.zeros((self.n, 1)) if μ_0 is None else μ_0
        self.Σ = np.eye(self.n) if Σ_0 is None else Σ_0

        # residuals 
        self.r = None

    def predict(self, t, u = np.array([0])):
        self.μ_pred = self.F[t] @ self.μ + self.B @ u
        self.Σ_pred = self.F[t] @ self.Σ @ self.F[t].T + self.Q

        return self.μ_pred, self.Σ_pred

    def update(self, z):
        
        self.r = z - self.H @ self.μ_pred # residual
        self.S = self.H @ self.Σ_pred @ self.H.T + self.R # residual covariance
        K = self.Σ_pred @ self.H.T @ np.linalg.inv(self.S) # Kalman gain matrix
        self.μ = self.μ_pred + K @ self.r # updated (a posteriori) state estimate
        I = np.eye(self.n) # n-size identity matrix
        self.Σ = (I - K @ self.H) @ self.Σ_pred # updated (a posteriori) estimate covariance
        	# (I - np.dot(K, self.H)).T) + np.dot(np.dot(K, self.R), K.T)
        
        return self.μ, self.Σ
    
    def log_likelihood(self):
        ll = self.r.T @ np.linalg.inv(self.S) @ self.r \
            + np.log(np.linalg.det(self.S)) \
            + self.m * np.log(2*np.pi)
        return np.sum( -0.5 * ll )
    
    def fullForward(self, measurements, u = np.array([0])):
        tables = []
        loglik = 0
        # (re)-initialization
        self.μ, self.Σ = self.μ_0, self.Σ_0
        for t,z in enumerate(measurements.T):
            self.predict(t,u=u)
            self.update(z)
            table=(self.μ_pred, self.Σ_pred, self.μ, self.Σ)
            tables.append(table)
            loglik += self.log_likelihood()
        return (tables, loglik)
    
    @staticmethod
    def Backward(tables,F):
        
        N = len(tables)
        μ_tT = [0]*N
        Σ_tT = [0]*N

        # initialization
        μ_tT[-1] = tables[-1][2]
        Σ_tT[-1] = tables[-1][3]
        
        # backward pass
        for t in range(N-2, -1, -1):
            μ_t = tables[t][2]
            Σ_t = tables[t][3]
            μ_t1t = tables[t+1][0]
            Σ_t1t = tables[t+1][1]
            
            J = Σ_t @ F[t].T @ np.linalg.inv(Σ_t1t)
            μ_tT[t] = μ_t + J @ (μ_tT[t+1] - μ_t1t) 
            Σ_tT[t] = Σ_t + J @ (Σ_tT[t+1] - Σ_t1t) @ J.T
        return (μ_tT, Σ_tT)