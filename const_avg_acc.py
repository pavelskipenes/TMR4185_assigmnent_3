
import numpy as np

def const_avg_acc(m, k, c, P, h, u0, udot0):
    # setter opp oppsett for om hele arrayet skal tas i en engang
     u_doubledot_0 = (P[0] - c*udot0 - u0*k)/m # fra equilibrium for t= 0 fra 3.79
     u1 = (P[1] + m * u_doubledot_0 + (4 * m / h + c) * udot0 + (4 * m / h**2 + 2 * c / h) * u0) / (4 * m / h**2 + 2 * c / h + k)#Hentet fra 3.80
     
     u = np.zeros(len(P)) #lage tomt u-array som kan endres ved indeksering
     u[0] = u0
     u[1] = u1

     #definere verdier som skal endres kontinuerlig i løkke
     u_doubledot_i = u_doubledot_0
     u_iplus1 = u1
     u_i = u0
     u_dot_i = udot0

     for i in range(len(P) - 1):
               
          u_iplus1 = (P[i + 1] + m * u_doubledot_i + (4 * m / h + c) * u_dot_i + (4 * m / h**2 + 2 * c / h) * u_i) / (4 * m / h**2 + 2 * c / h + k)
          
          u[i + 1] = u_iplus1

          u_doubledot_iplus1 = 4/h**2 * (u_iplus1 - u_i - h*u_dot_i) - u_doubledot_i
               
          u_dot_iplus1 = 2/h * (u_iplus1 - u_i) - u_dot_i

          # definerer nye i verdier basert på i+1 fra forrige i
          u_doubledot_i = u_doubledot_iplus1
          u_i = u_iplus1
          u_dot_i = u_dot_iplus1
     
     return u