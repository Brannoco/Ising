
# coding: utf-8

# In[56]:


# Code Mark 1.7: Ising Model in 2D

# Importing all the necessary packages I will use:
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import numpy.random as rnd

# Before I do anything I need to make sure the specifics of the size of the lattice as a 12x12 sized, number of steps, number of temperature points used, number of points in order to reach an equilibrium state are all accounted for
num=240
steps_monte=150
temp_points=100
steps_equil=150
num_a=1/(steps_monte*(num**2))
num_b=1/((steps_monte**2)*(num**2))

# Next I must create an intitial spin state for my molecule: with 2 different values e.g. (-1,1) or (0,1) etc.
def startspin(num):   
    spin=rnd.randint(2,size=(num,num))-1
    return spin

# Now I need to define the energy function. This will give the energy of any given spin state:
def Energy(Q):
    starting_energy=0
    for i in range(len(Q)):
        for j in range(len(Q)):
            g=Q[i,j]
            n_y=Q[(i+1)%num,j]+Q[i,(j+1)%num]+Q[(i-1)%num,j]+Q[i,(j-1)%num]
            starting_energy+=g*-n_y
    return starting_energy/4

# Now I need to code the Monte Carlo Steps with Metro ALgorithm for the ising model:
def montestep(Q,P):
    for i in range(num):
        for j in range(num):
                x=rnd.randint(0,num)
                y=rnd.randint(0,num)
                z=Q[x,y]
                n_y=Q[(x+1)%num,y]+Q[x,(y+1)%num]+Q[(x-1)%num,y]+Q[x,(y-1)%num]
                l=2*z*n_y
                if l<0:
                    z*=-1
                elif rnd.rand()<np.exp(P*-l):
                    z*=-1
                Q[x,y]=z
    return Q

#Now I need to define the magnetization function. This gives the Magnetization of any given spin state:
def Magnetization(Q):
    magnetization=np.sum(Q)
    return magnetization

# Here I am saying that the temperature points need to have a midpoint where we expect the disorder to begin and to give us random points between this value as 1<T<5
temp_mid = 2.25;    
Temp=rnd.normal(temp_mid,0.5,temp_points)
Temp=Temp[(Temp>0.5)&(Temp<5)]   
temp_points=np.size(Temp)

# I need to define 0 points for each of the 4 physical quantities I wish to find:
Magnetization=np.zeros(temp_points)
SpecificHeat=np.zeros(temp_points)  
Energy=np.zeros(temp_points)
Susceptibility=np.zeros(temp_points)

# Now I need to implement each function above by creating an Ising Function with them:
for i in range(len(Temp)):
    E_a=M_a=0
    E_b=M_b=0
    Q=startspin(num)
    init_Temp_a=1/Temp[i] 
    init_Temp_b=init_Temp_a**2

    for j in range(steps_monte):
        montestep(Q,init_Temp_a)            
        Energy_calc=Energy[Q]  
        Mag_calc=Magnetization[Q]       
        E_a=E_a+(Energy_calc)
        M_a=M_a+(Mag_calc)
        M_b=M_b+(Mag_calc**2)
        E_b=E_b+(Energy_calc**2)
        Magnetization[j]=num_a*M_a[i][j]
        SpecificHeat[j]=(num_a*E_b[i][j]*(E_a[i][j]**2))*init_Temp_b
        Susceptibility[j]=(num_a*M_b[i][j]*(M_a[i][j]**2))*init_Temp_a
        Energy[j]=num_a*E_a[i][j]

    for k in range(steps_equil):        
        montestep(Q,init_Temp_a)   
    
#Finally I can plot the 4 graphs I need to show the monte carlo steps of the ising model and can then conduct my research thereafter.
        
plt.plot(Temp, Energy,'d', color="#8A2BE2")
plt.xlabel("Temperature (Kelvin)")
plt.ylabel("Energy (Arbitrary Units)")
plt.show()

plt.plot(Temp, abs(Magnetization),'x', color="#7FFF00")
plt.xlabel("Temperature (Kelvin)")
plt.ylabel("Magnetization (Arbitrary Units)")
plt.show()

plt.plot(Temp, SpecificHeat, 'd', color="#00FFFF")
plt.xlabel("Temperature (Kelvin)")
plt.ylabel("Specific Heat (Arbitrary Units)")
plt.show()

plt.plot(Temp, Susceptibility, 'x', color="#0000FF")
plt.xlabel("Temperature (Kelvin)")
plt.ylabel("Susceptibility (Arbitrary Units)")
plt.show()










