import numpy as np

def state_space(m,k,c,P,h,u0,udot0):
    n=len(P)    #fordi P er en n lang vektor
    
    u=np.zeros(n)
    udot=np.zeros(n)

    u[0]=u0
    udot[0]=udot0
    
    x=np.array([u,udot])
    A=np.array([[0,1],[-k/m,-c/m]])   #matriseform for newton 2.L
    B=np.array([[0],[1]])               #bare mhp last for x_2 derivert->akselerasjon

    #tar ikke hensyn til C, ettersom forholdstallet bare er 1.

    for i in range(1, n):
        dx=A@x + B*P[i-1]
        x=x+dx*h

        u[i]=x[0][0]
        udot[i]=x[0][1]

    return u



        
    
