import numpy as np
import matplotlib.pyplot as plt


##########################
#w_0=4*w_p, undamped

def u(t,k,P_0,w_p):
    return (P_0 / k) * (1/2 
        + (1 / (np.pi * 15)) * np.sin(w_p * t) 
        + (1 / (3 * np.pi * 7)) * np.sin(3 * w_p * t) 
        - (1 / (5 * np.pi * 9)) * np.sin(5 * w_p* t) 
        - (1 / (7 * np.pi * 33)) * np.sin(7 * w_p * t)
    )

    #Hvordan skrive som sum av n ganger?

time=np.linspace(0,51,300)

u_t=[]
u_t_ud=[]

k=10   #N/m
P_0=10  #N
w_p=2*np.pi/10  #rad/sec

m=2 #kg
w_0=np.sqrt(k/m)    #rad/sec
w_p_new=w_0/4

for t in time:
    u_t.append(u(t,k,P_0,w_p_new))

def beta_n(n,w_p):
    return n*w_p/w_0

def u_ud(t,k,P_0,w_p,n):
    result=P_0/(2*k)

   
    for i in range(0,n):
        if i%2==0:
            result+=0

        else:
            result+=(1/k)*1/(1-beta_n(i,w_p)**2)*2*P_0/(i*np.pi)
    return result

""" return_value = sum([const_factor * i for i in range(0, 10) if i % 2 == 1]) """
n=1000

for t in time:
    u_t_ud.append(u_ud(t,k,P_0,w_p,n))

P_load=[]

def P(t, P0, w_p, n_terms):
    res = P0/2
    for n in range(1, n_terms + 1, 2):
        res += P0 * (2 / (n * np.pi)) * np.sin(n * w_p * t)
    
    return res

for i in time:
    P_load.append(P(i,P_0,w_p,4))

plt.plot(time,u_t)
plt.xlabel("Time [s]")
plt.ylabel("Displacement [m]")
plt.show()

plt.plot(time,P_load)
plt.xlabel("Time [s]")
plt.ylabel("Force [N]")
plt.show()

k_temp=1
c_temp=0.1
delta_u=1

k_list=np.linspace(k_temp,20)
c_list=np.linspace(c_temp,20)

k_values=[]
c_values=[]
k_res=15.53
c_res=10.69

def kc_check():
    min_k=0
    min_c=0
    for i in k_list:
        temp=(P_0/i)/(np.sqrt((1-(2*np.pi/(10*np.sqrt(i/m)))**2)**2 + (2 *c_temp/(2*np.sqrt(m*i)) * (2*np.pi/(10*np.sqrt(i/m))))**2))
        k_values.append(temp)
        if temp>=delta_u:
            min_k=i
            k_res=i

    for j in c_list:
        temp=(P_0/k_temp)/(np.sqrt((1-(2*np.pi/(10*np.sqrt(k_temp/m)))**2)**2 + (2 *j/(2*np.sqrt(m*k_temp)) * (2*np.pi/(10*np.sqrt(k_temp/m))))**2))
        c_values.append(temp)
        if temp>=delta_u:
            min_c=j
            c_res=j
    
    return (f'Min k when c=0.1: {min_k:.2f} [N/m], Min c when k=1: {min_c:.2f} [Nm/s]')

print(kc_check())

plt.plot(k_list,k_values, label="A(k), c=0.1 [Nm/s]")
plt.plot(c_list,c_values, label="A(c), k=1 [N/m]")
plt.scatter(k_res, 1, color='blue')
plt.plot(c_res, 1, 'ro')
plt.text(19, 4, f' K_min = {k_res:.2f} ', fontsize=12, color='blue', ha='right')
plt.text(12, 4, f' C_min = {c_res:.2f} ', fontsize=12, color='red', ha='right')
plt.axhline(y=1, color='green', linestyle='--', label="Amplitude = 1 [m]")
plt.legend(loc='upper right')
plt.xlabel("K [N/m] , C [Nm/s]")
plt.ylabel("Amplitude [m]")
plt.show()


