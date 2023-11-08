import numpy as np

class KalmanFilter(object):
    def __init__(self, F = None, B = None, H = None, Q = None, R = None, R_a = None, R_b = None, R_c = None, μ_0 = None, Σ_0 = None):

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

    def predict(self, t=0, u = np.array([0])):
        if len(self.F.shape) == 3:
            self.μ_pred = self.F[t] @ self.μ + self.B @ u
            self.Σ_pred = self.F[t] @ self.Σ @ self.F[t].T + self.Q
        else:
            self.μ_pred = self.F @ self.μ + self.B @ u
            self.Σ_pred = self.F @ self.Σ @ self.F.T + self.Q

        return self.μ_pred, self.Σ_pred

    def update(self, z):

        if np.isnan(z).all():
            self.μ = self.μ_pred
            self.Σ = self.Σ_pred
            
        else:
            self.r = z - self.H @ self.μ_pred # residual
            self.S = self.H @ self.Σ_pred @ self.H.T + self.R # residual covariance
            K = self.Σ_pred @ self.H.T @ np.linalg.inv(self.S) # Kalman gain matrix
            self.μ = self.μ_pred + K @ self.r # updated (a posteriori) state estimate
            I = np.eye(self.n) # n-size identity matrix
            self.Σ = (I - K @ self.H) @ self.Σ_pred # updated (a posteriori) estimate covariance
        
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
            
            if len(F.shape) == 3:
                J = Σ_t @ F[t].T @ np.linalg.inv(Σ_t1t)
            else:
                J = Σ_t @ F.T @ np.linalg.inv(Σ_t1t)
            μ_tT[t] = μ_t + J @ (μ_tT[t+1] - μ_t1t) 
            Σ_tT[t] = Σ_t + J @ (Σ_tT[t+1] - Σ_t1t) @ J.T
        return (μ_tT, Σ_tT)
    


class KalmanFilter_rotating_Q(object):
    def __init__(self, F = None, B = None, H = None, Q = None, R = None, R_a = None, R_b = None, R_c = None, μ_0 = None, Σ_0 = None):

        if(F is None or H is None):
            raise ValueError("Set proper system dynamics.")
        
        # Parameters of the model
        self.m, self.n = H.shape # measurement and hidden state dimensions
        self.F = F # transition matrices
        self.H = H # Observation matrix
        self.B = np.zeros(1) if B is None else B # control matrix (optional)
        self.Q = np.eye(self.n) if Q is None else Q # transition noise covariance
        self.σ_r = self.Q[0,0]
        self.σ_φ = self.Q[1,1]
        self.R = np.eye(self.m) if R is None else R # observation noise covariance
        self.S = None

        # exponential decay of R as a function of z :  R(z) = a * exp(-b * z) + c
        self.R_a = 0.01 if R_a is None else R_a
        self.R_b = 0.7  if R_b is None else R_b
        self.R_c = 0.2  if R_c is None else R_c

        # save initial state for backward pass
        self.μ_0 = μ_0
        self.Σ_0 = Σ_0
        
        # predicted state
        self.μ_pred = None
        self.Σ_pred = None        

        # updated / current state
        self.μ = np.zeros((self.n, 1)) if μ_0 is None else μ_0
        self.Σ = np.eye(self.n) if Σ_0 is None else Σ_0
        self.r = np.sqrt( self.μ[0] ** 2 + self.μ[1] ** 2)
        self.φ = np.arctan2(self.μ[1], self.μ[0])

        # residuals 
        self.residual = None
    
    def set_polar_coordinates(self,μ):
        self.r = np.sqrt( μ[0] ** 2 + μ[1] ** 2)
        self.φ = np.arctan2(μ[1], μ[0])

    def setQ(self):
        self.Q[0,0] = self.σ_r * np.cos(self.φ) ** 2 + self.σ_φ * np.sin(self.φ) ** 2
        self.Q[0,1] = (self.σ_r - self.σ_φ) * np.sin(self.φ) * np.cos(self.φ)
        self.Q[1,0] = self.Q[0,1]
        self.Q[1,1] = self.σ_r * np.sin(self.φ) ** 2 + self.σ_φ * np.cos(self.φ) ** 2

    def setR(self, z):
        self.R = np.diag( self.R_a * np.exp(-self.R_b * z) + self.R_c )

    def predict(self, t=0, u = np.array([0])):
        self.μ_pred = self.F @ self.μ + self.B @ u
        self.set_polar_coordinates(self.μ_pred)
        self.setQ()
        self.Σ_pred = self.F @ self.Σ @ self.F.T + self.Q

        return self.μ_pred, self.Σ_pred

    def update(self, z):

        if np.isnan(z).all():
            self.μ = self.μ_pred
            self.Σ = self.Σ_pred
            
        else:
            self.residual = z - self.H @ self.μ_pred # residual
            self.setR(z)
            self.S = self.H @ self.Σ_pred @ self.H.T + self.R # residual covariance
            K = self.Σ_pred @ self.H.T @ np.linalg.inv(self.S) # Kalman gain matrix
            self.μ = self.μ_pred + K @ self.residual # updated (a posteriori) state estimate
            I = np.eye(self.n) # n-size identity matrix
            self.Σ = (I - K @ self.H) @ self.Σ_pred # updated (a posteriori) estimate covariance
        
        return self.μ, self.Σ
    
    def log_likelihood(self):
        ll = self.residual.T @ np.linalg.inv(self.S) @ self.residual \
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
            
            if len(F.shape) == 3:
                J = Σ_t @ F[t].T @ np.linalg.inv(Σ_t1t)
            else:
                J = Σ_t @ F.T @ np.linalg.inv(Σ_t1t)
            μ_tT[t] = μ_t + J @ (μ_tT[t+1] - μ_t1t) 
            Σ_tT[t] = Σ_t + J @ (Σ_tT[t+1] - Σ_t1t) @ J.T
        return (μ_tT, Σ_tT)
    

class KalmanFilter_k_hidden(object):
    def __init__(self, F = None, B = None, H = None, Q = None, R_a = None, R_b = None, R_c = None, μ_0 = None, Σ_0 = None, dx=None):

        if(F is None or H is None):
            raise ValueError("Set proper system dynamics.")
        
        # Parameters of the model
        self.m, self.n = H.shape # measurement and hidden state dimensions
        self.F = F # transition matrices
        self.H = H # Observation matrix
        self.B = np.zeros(1) if B is None else B # control matrix (optional)
        self.Q = np.eye(self.n) if Q is None else Q # transition noise covariance
        self.σ_r = self.Q[0,0]
        self.σ_φ = self.Q[1,1]
        self.R = np.eye(self.m) # observation noise covariance
        self.S = None

        # exponential decay of R as a function of z :  R(z) = a * exp(-b * z) + c
        self.R_a = 0.008 if R_a is None else R_a
        self.R_b = 0.8  if R_b is None else R_b
        self.R_c = 0.05  if R_c is None else R_c

        # save initial state for backward pass
        self.μ_0 = μ_0
        self.Σ_0 = Σ_0
        
        # predicted state
        self.μ_pred = None
        self.Σ_pred = None

        # updated / current state
        self.μ = np.zeros((self.n, 1)) if μ_0 is None else μ_0 # current state (a, b, a', b')
        self.Σ = np.eye(self.n) if Σ_0 is None else Σ_0 # current covariance matrix
        self.r = np.sqrt( self.μ[0] ** 2 + self.μ[1] ** 2) # current radius [bp]
        self.φ = np.arctan2(self.μ[1], self.μ[0]) # current angle [rad]
        self.dφ = (self.μ[3]*self.μ[0] - self.μ[2]*self.μ[1])/self.r**2 # angular velocity [rad/bp]
        self.dx = dx # distance between bins [bp]
        self.θ = self.dφ * self.dx # rotation angle [rad]

        # residuals
        self.residual = None
    
    def set_polar_coordinates(self,μ):
        self.r = np.sqrt( μ[0] ** 2 + μ[1] ** 2)
        self.φ = np.arctan2(μ[1], μ[0])
        self.dφ = (μ[3]*μ[0] - μ[2]*μ[1])/self.r**2
        self.θ = self.dφ * self.dx

    def setF(self):
        self.F[0,0] = np.cos(self.θ)
        self.F[0,1] = -np.sin(self.θ)
        self.F[1,0] = np.sin(self.θ)
        self.F[1,1] = np.cos(self.θ)

    def setQ(self):
        self.σ_φ = (self.θ/10)**2
        self.Q[0,0] = self.σ_r * np.cos(self.φ) ** 2 + self.σ_φ * np.sin(self.φ) ** 2
        self.Q[0,1] = (self.σ_r - self.σ_φ) * np.sin(self.φ) * np.cos(self.φ)
        self.Q[1,0] = self.Q[0,1]
        self.Q[1,1] = self.σ_r * np.sin(self.φ) ** 2 + self.σ_φ * np.cos(self.φ) ** 2

    def setR(self, z):
        self.R = np.diag( self.R_a * np.exp(-self.R_b * z) + self.R_c )

    def predict(self, t=0, u = np.array([0])):
        self.setF()
        self.μ_pred = self.F @ self.μ + self.B @ u
        self.set_polar_coordinates(self.μ_pred)
        self.setQ()
        self.Σ_pred = self.F @ self.Σ @ self.F.T + self.Q

        return self.μ_pred, self.Σ_pred

    def update(self, z):

        if np.isnan(z).all():
            self.μ = self.μ_pred
            self.Σ = self.Σ_pred
            
        else:
            self.residual = z - self.H @ self.μ_pred # residual
            self.setR(z)
            self.S = self.H @ self.Σ_pred @ self.H.T + self.R # residual covariance
            K = self.Σ_pred @ self.H.T @ np.linalg.inv(self.S) # Kalman gain matrix
            self.μ = self.μ_pred + K @ self.residual # updated (a posteriori) state estimate
            I = np.eye(self.n) # n-size identity matrix
            self.Σ = (I - K @ self.H) @ self.Σ_pred # updated (a posteriori) estimate covariance
        
        return self.μ, self.Σ
    
    def log_likelihood(self):
        ll = self.residual.T @ np.linalg.inv(self.S) @ self.residual \
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
            
            if len(F.shape) == 3:
                J = Σ_t @ F[t].T @ np.linalg.inv(Σ_t1t)
            else:
                J = Σ_t @ F.T @ np.linalg.inv(Σ_t1t)
            μ_tT[t] = μ_t + J @ (μ_tT[t+1] - μ_t1t) 
            Σ_tT[t] = Σ_t + J @ (Σ_tT[t+1] - Σ_t1t) @ J.T
        return (μ_tT, Σ_tT)