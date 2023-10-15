import numpy as np

def initial_conditions(n, L, T, dt):
    # number of particles
    # size of the box (and frame)

    index = np.linspace(1,n,n)
    # mass of argon
    mass = np.ones(n) * 6.6335209e-26
    
#     T is the sum of the energies of the particles
    speed = np.cbrt(np.sqrt((2*T)/(6.6335209e-26*n)))
    

    pos = np.random.uniform(-L/2, L/2, n*3)
    vel = ((-1)**np.random.randint(0,2,n*3)) * np.ones(n*3) * speed


    pos = [0.5e-5,0,0,-0.5e-5,0,0]
    vel = np.zeros(6)
    file = open(f"cache/initial_conditions_{n}n_{L}m_{dt}s.csv", "w+")

    for i in range(0,n):
        string = ["Particle "+str(int(index[i]))+"\n",str(mass[i])+"\n",str(pos[3*i])+", "+str(pos[3*i+1])+", "+str(pos[3*i+2])+"\n",str(vel[3*i])+", "+str(vel[3*i+1])+", "+str(vel[3*i+2])+"\n"]
        file.writelines(string)
    file.close()
