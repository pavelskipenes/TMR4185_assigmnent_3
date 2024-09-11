#rasnlch
import numpy as np
def const_avg_acc(m, k , c ,P, h , u0 , udot0):
    u_doubledot_0 = ((P[h] - c * udot0 - k * u0)/m) # fra equilibrium for t= 0 fra 3.79
    u = ((P[h + 1] + m*u_doubledot_0 
         + (4*m/h + c)*udot0 
         + (4*m/h**2 + 2*c/h)*u0)
         /(4*m/h**2 + 2*c/h + k))#Hentet fra 3.80
    return u


def ny_const(m, k, c, P, h, u0, udot0):
    # setter opp oppsett for om hele arrayet skal tas i en engang
     u_doubledot_0 = ((P[0] - c * udot0 - k * u0)/m) # fra equilibrium for t= 0 fra 3.79
     u1 = ((P[1] + m*u_doubledot_0 
          + (4*m/h + c)*udot0 
          + (4*m/h**2 + 2*c/h)*u0)
          /(4*m/h**2 + 2*c/h + k))#Hentet fra 3.80
     u = np.zeros(len(P)) #lage tomt u-array som kan endres ved indeksering


     for i in range(len(P) - 1):#vil ikke indeksere utenfor rammer med i + 1
          if i == 1:#Kan være dette må være 0 og neste større enn null, 
                    #men forstår ikke hvordan jeg skal få verdier for 0
               u_doubledot_i = u_doubledot_0
               u_iplus1 = u1
               u_i = 0
               u_dot_i = udot0
               #Når vi finner rett grense kan dette settes før løkke og heller endre range (fjerne if)

          elif i > 1: # definerer nye i verdier basert på i+1 fra forrige i
               u_doubledot_i = u_doubledot_iplus1
               u_i = u_iplus1
               u_dot_i = u_dot_iplus1
               
               #kjører utregning i henhold til metode
               u_doubledot_iplus1 = (2*h*(P[i+1] - k*u_iplus1 - c*u_i - c*u_doubledot_i / 2*h)
                                /(2*h*m + c)) #Mishandling fra å sette 3.78 inn i 3.79
               u_dot_iplus1 = u_dot_i + 1/2 * (u_doubledot_i + u_doubledot_iplus1)*h
               u_iplus1 = ((P[i+1] + m*u_doubledot_i 
                            + (4*m/h + c)*u_dot_i 
                            + (4*m/h**2 + 2*c/h)*u_i)
                            /(4*m/h**2 + 2*c/h + k))
          u[i + 1] = u_iplus1
     
     return u