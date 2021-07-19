def acoustic_impedance(v,rho):
    return v*rho

def poisson_ratio(vp,vs):
    return (vp**2 - 2*(vs**2))/(2*(vp**2 - vs**2))
    
class Zoeppritz_approximation:

    def __init__(self,theta,vp_1,vp_2,vs_1,vs_2,rho_1,rho_2):
        self.theta = theta
        self.vp_1 = vp_1
        self.vp_2 = vp_2
        self.vs_1 = vs_1
        self.vs_2 = vs_2
        self.rho_1 = rho_1
        self.rho_2 = rho_2
        self.delta_vp = self.vp_2-self.vp_1
        self.delta_vs = self.vs_2-self.vs_1
        self.delta_rho = self.rho_2-self.rho_1
        self.vp = (self.vp_1+self.vp_2)/2
        self.vs = (self.vs_1+self.vs_2)/2
        self.rho = (self.rho_1+self.rho_2)/2
        
        self.AI_1 = acoustic_impedance(self.vp_1,self.rho_1)
        self.AI_2 = acoustic_impedance(self.vp_2,self.rho_2)
        
        self.poisson_ratio_1 = poisson_ratio(self.vp_1,self.vs_1)
        self.poisson_ratio_2 = poisson_ratio(self.vp_2,self.vs_2)
        
        self.A = (1/2)*((self.delta_vp/self.vp)+(self.delta_rho/self.rho))
        self.B = (self.delta_vp/(2*self.vp)) - 4*((self.vs/self.vp)**2)*(self.delta_vs/self.vs) - 2*((self.vs/self.vp)**2)*(self.delta_rho/self.rho)
        self.C = (1/2)*(self.delta_vp/self.vp)
        
    def aki_and_richards(self):
        return self.A + self.B*(np.sin(self.theta)**2) + self.C*((np.sin(self.theta)**2))*((np.tan(self.theta)**2))
    
    def shuey(self):
        return self.A + self.B*(np.sin(self.theta)**2)
    
    def hilterman(self):
        mean_poisson = (self.poisson_ratio_2+self.poisson_ratio_1)/2
        return ((self.AI_2 - self.AI_1)/(self.AI_2+self.AI_1))*(np.cos(self.theta)**2) + ((self.poisson_ratio_2-self.poisson_ratio_1)/((1-mean_poisson)**2))*(np.sin(self.theta)**2)