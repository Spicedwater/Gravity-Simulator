import numpy as np
import math
from matplotlib import pyplot as plt
from matplotlib import animation


D = 2
G = 6.67428e-11
N = 2
AU = 1.496e11
scale = 200/AU

class Body():
    def __init__(self, name, mass, p, v, colour):
        self.name = name
        self.mass = mass
        self.p = p
        self.v = v
        self.F = np.array([])
        self.colour = colour

    def force(self, other):
        if self is not other:
            dx = other.p[0] - self.p[0]
            dy = other.p[1] - self.p[1]

            r = (dx**2 + dy**2) ** 0.5

            if r == 0:
                raise ValueError("Collision between objects {} and {}"
                                 .format(self.name, other.name))

            f = G * self.mass * other.mass / r**2

            ang = math.atan2(dy, dx)
            fx = f * math.cos(ang)
            fy = f * math.sin(ang)

            return [fx, fy]

    def tot_force(self, bodies):
        fx_tot = fy_tot = 0
        for body in bodies:
                if body is self:
                    continue
                force = self.force(body)
                fx_tot += force[0]
                fy_tot += force[1]
            
        self.F = [fx_tot, fy_tot]


class Rocket(Body):
    def __init__(self, name, mass, origin, target):
        self.name = name
        self.mass = mass
        self.origin = origin
        self.target = target

    def hohmann(self, bodies):
        r = (self.origin.p[0]**2 + self.origin.p[1]**2) ** 0.5 + (self.target.p[0]**2 + self.target.p[1]**2) ** 0.5 / 2

        self.tot_force(bodies)

        T = (4 * math.pi**2 * self.mass * r / self.F) ** 0.5

        mass_sum = 0
        for body in bodies:
            mass_sum += body.mass

        r1 = (self.origin.p[0]**2 + self.origin.p[1]**2)**0.5
        r2 = (self.target.p[0]**2 + self.target.p[1]**2)**0.5

        dv1 = (G * mass_sum / r1)**0.5 * ((2 * r2 / (r1 + r2))**0.5 - 1)
        dv2 = (G * mass_sum / r2)**0.5 * (1 - (2 * r1 / (r1 + r2))**0.5)

        



def print_data(step, bodies):
    print("Day #{}".format(step))
    for body in bodies:
        data = '{:<8}  Pos.={:>6.2f} {:>6.2f} Vel.={:>10.3f} {:>10.3f}'.format(
           body.name, body.p[0]/AU, body.p[1]/AU, body.v[0]/1000, body.v[1]/1000)
        print(data)
    print()


def refresh(bodies, dt):
    #dt = 24 * 3600
    #elapsed_time = 0

    #while elapsed_time < (100 * dt):
        #elapsed_time += dt
    #print_data((elapsed_time//dt), bodies)
    
    for body in bodies:
        body.tot_force(bodies)

    for body in bodies:
        for i in range(0, 2):
            #print(body.v[i])
            body.v[i] += body.F[i] / body.mass * dt
            #print(body.F[i], body.name)
            #print(body.F[i] / body.mass * dt)
            #print(body.v[i])
            body.p[i] += body.v[i] * dt



bodies = []

# for x in range(0, N):
#     print("\nEnter the mass of body {} in kg:".format(x + 1))
#     mass = float(input(">>>   "))

#     print("Enter the p of body {} in km in the format [x, y]".format(x + 1))
#     p = np.array(input(">>>   "))

#     print("Enter the v of body {} in m/s in the format [x, y]".format(x + 1))
#     v = np.array(input(">>>   "))

    # print("Enter the a of body {} in km in the format [x, y, z]".format(x + 1))
    # a = np.array(input(">>>   "))


# sun = Body('Sun', 1.989e30, [0, 0], [0, 0], 'yellow')
# bodies.append(sun)

# mercury = Body('Mercury', 3.33011e23, [6.982e10, 0], [0, -38860], 'orange')
# bodies.append(mercury)

# venus = Body('Venus', 4.8675e24, [1.0894e11, 0], [0, -34790], 'grey')
# bodies.append(venus)

# earth = Body('Earth', 5.972e24, [1.522e11, 0], [0, -29290], 'blue')
# bodies.append(earth)

earth2 = Body('Earth2', 5.972e24, [0, 0], [-1, -11.7], 'blue')
bodies.append(earth2)

moon = Body('Moon', 7.346e22, [earth2.p[0] + 3.544e8, 0], [0, 970], 'DarkGrey')
bodies.append(moon)

# mars = Body('Mars', 6.4171e23, [2.4923e11, 0], [0, -21970], 'red')
# bodies.append(mars)

# jupiter = Body('Jupiter', 1.89819e27, [8.1662e11, 0], [0, -12440], 'orange')
# bodies.append(jupiter)

#refresh(bodies)

fig = plt.figure()
#plt.axis((-6, 6, -6, 6))

dt = 4 * 3600
elapsed_time = 0
while True:
    elapsed_time += dt
    print_data(elapsed_time/(24 * 3600), bodies)
    #plt.clf()
    refresh(bodies, dt)
    for body in bodies:
        plt.plot(body.p[0]/AU, body.p[1]/AU, color=body.colour, marker='.', ms=2)
        plt.pause(0.0001)
    # plt.plot(earth.p[0]/AU, earth.p[1]/AU, color=earth.colour, marker='.', ms=2)
    # plt.plot(earth2.p[0]/AU, earth2.p[1]/AU, color=earth2.colour, marker='.', ms=2)
    # plt.plot([earth.p[0]/AU, earth2.p[0]/AU], [earth.p[1]/AU, earth2.p[1]/AU], linewidth=0.5, color='blue')
    # plt.plot(sun.p[0]/AU, sun.p[1]/AU, color=sun.colour, marker='.', ms=2)
    # plt.pause(0.1)


        


