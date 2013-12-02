import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def run_md(steps = 1000000, dt = 0.002, n = 50, sigma = 1.0):
    SAVE_INTERVAL = 10000000
    rc = 1
    #set positions and velocities
    particles = np.random.uniform(0, 9*sigma, (n,3))
    velocities = np.random.uniform(-1, 1, (n,3))
    #ensure center of mass has zero velocity
    for i in range(3):
        excess_v = np.sum(velocities[:,i])
        velocities[:,i] = velocities[:,i]-excess_v/n
    fig = plt.figure()
    fig2 = plt.figure()
    axis = fig.add_subplot(111, projection='3d')
    axis2 = fig2.add_subplot(111)
    #calc initial accel
    accelerations = np.zeros((n,3))
    expanded_particles = np.zeros((4*n,3))
    for i in range(n):
        particle = particles[i,:]
        expanded_particles[4*i,:] = particle
        for j in range(3):
            unit_v = np.zeros((1,3))
            unit_v[0,j] = 1
            if particle[j] > 9*sigma/2:
                expanded_particles[4*i+j+1,:] = particle - 9*sigma*unit_v
            else:
                expanded_particles[4*i+j+1,:] = particle + 9*sigma*unit_v
    #2nd: find each particles nearest neighbors:
    nearest_neighbors = []
    for i in range(n):
        nearest_neighbors.append([])
        p = particles[i,:]
        for j in range(4*n):
            test_p = expanded_particles[j,:]
            dx = p[0] - test_p[0]
            dy = p[1] - test_p[1]
            dz = p[2] - test_p[2]
            if (np.power(dx*dx + dy*dy + dz*dz, 0.5) < rc) & (j != 4*i):
                nearest_neighbors[-1].append(test_p)
    #3rd calc resultant acceleration per particle:
    for i in range(n):
        particle = particles[i,:]
        ax = 0.0
        ay = 0.0
        az = 0.0
        for neighbor in nearest_neighbors[i]:
            x = particle[0] - neighbor[0]
            y = particle[1] - neighbor[1]
            z = particle[2] - neighbor[2]
            r = np.power(x*x + y*y + z*z, 0.5)
            ax = ax + x*(1-1/r)
            ay = ay + y*(1-1/r)
            az = az + z*(1-1/r)
        accelerations[i,:] = np.array([ax, ay, az])
    potential_energies = np.zeros((steps,1))
    kinetic_energies = np.zeros((steps,1))
    for step in range(steps):
        particles = particles + dt*velocities + dt*dt*accelerations/2
        #check if out of bounds
        for i in range(n):
            for j in range(3):
                if particles[i,j] < 0:
                    particles[i,j] = 9*sigma+particles[i,j];
                elif particles[i,j] > 9*sigma :
                    particles[i,j] = 9*sigma - particles[i,j];
        #calc accelerationnn
        #1st: expand box
        expanded_particles = np.zeros((4*n,3))
        for i in range(n):
            particle = particles[i,:]
            expanded_particles[4*i,:] = particle
            for j in range(3):
                unit_v = np.zeros((1,3))
                unit_v[0,j] = 1
                if particle[j] > 9*sigma/2:
                    expanded_particles[4*i+j+1,:] = particle - 9*sigma*unit_v
                else :
                    expanded_particles[4*i+j+1,:] = particle + 9*sigma*unit_v
        #2nd: find each particles nearest neighbors:
        nearest_neighbors = []
        for i in range(n):
            nearest_neighbors.append([])
            p = particles[i,:]
            for j in range(4*n):
                test_p = expanded_particles[j,:]
                dx = p[0] - test_p[0]
                dy = p[1] - test_p[1]
                dz = p[2] - test_p[2]
                if (np.power(dx*dx + dy*dy + dz*dz, 0.5) < rc) & (j != 4*i):
                    nearest_neighbors[-1].append(test_p)
        #3rd calc resultant acceleration per particle:
        new_accel = np.zeros((n,3))
        potential_energy = 0
        for i in range(n):
            ax = 0.0
            ay = 0.0
            az = 0.0
            particle = particles[i,:]
            for neighbor in nearest_neighbors[i]:
                x = particle[0] - neighbor[0]
                y = particle[1] - neighbor[1]
                z = particle[2] - neighbor[2]
                r = np.power(x*x + y*y + z*z, 0.5)
                ax = ax + x*(1-1/r)
                ay = ay + y*(1-1/r)
                az = az + z*(1-1/r)
                potential_energy = potential_energy + pow(1-r, 2)/2
            new_accel[i,:] = np.array([ax, ay, az])
        velocities = velocities + (new_accel+accelerations)*dt/2
        accelerations = new_accel
        #calc total energy
        kinetic_energy = 0
        for i in range(n):
            kinetic_energy = kinetic_energy + np.sum(velocities[i,:]*velocities[i,:])/2
        #kinetic_energy = np.sum([np.sum(velocities[i,:]*velocities[i,:]) for i in range(n)])/2
        potential_energies[step] = potential_energy
        kinetic_energies[step] = kinetic_energy
        #plot
        if step % SAVE_INTERVAL == 0:
            axis.clear()
            axis.set_xlim((0,9*sigma))
            axis.set_ylim((0,9*sigma))
            axis.set_zlim((0,9*sigma))
            axis.grid(False)
            axis.scatter(particles[:,0], particles[:,1], particles[:,2])
            plt.savefig('./figs/fig_' + str(step / SAVE_INTERVAL) + '.png')
        print step
    axis2.plot(np.linspace(1, steps, steps), kinetic_energies, 'b-', label='kinetic energy')
    axis2.plot(np.linspace(1, steps, steps), potential_energies, 'g-', label='potential energy')
    axis2.plot(np.linspace(1, steps, steps), potential_energies + kinetic_energies, 'r-', label='total energy')
    axis2.set_xlabel('time')
    axis2.set_ylabel('energy')
    axis2.legend(loc=7)
    plt.savefig('./figs/energies.png')

if __name__ == "__main__":
    run_md()
