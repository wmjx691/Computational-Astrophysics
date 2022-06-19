import numpy as np

class PhysicsBody:

    def __init__(self, initPosx, initPosy,
                 initAngle, initVelocity):
        """
        Initialize the initial position, angle, and velocity
        """
        self.mass=1.0
        self.reset(initPosx, initPosy,
                   initAngle, initVelocity)

        return

    def reset(self, initPosx, initPosy,
              initAngle, initVelocity):
        """
        Reset the postion and angle for replay
        """
        angle = np.pi * initAngle / 180.0
        self.posx = initPosx
        self.posy = initPosy
        self.velx = initVelocity * np.cos(angle)
        self.vely = initVelocity * np.sin(angle)
        return

    def add_gravity(self, gravity = -980.0):
        """
        setup the gravity
        """
        self.gravity = gravity
        return

    def get_acceleration(self):
        """
        Get the acceleration

        output: the acceleration in x- and y- directions
                ax [float]
                ay [float]
        """

        # gravitational forces


        # other forces


        # calculate the acceleration
        
        return ax, ay

    def update(self, dt=0.01):
        """
        Do one update with dt = 0.01 [sec]

        """
        posx = self.posx
        posy = self.posy
        velx = self.velx
        vely = self.vely

        yin = np.array([posx,
                        posy,
                        velx,
                        vely])

        def derive(yin):
            """
            Input:
                yin : numpy array

            Output:
                yout : numpy array : the dydt
            """
            dydt0 = yin[2]      # dxdt = vx
            dydt1 = yin[3]      # dydt = vy
            dydt2 = 0      # dvxdt = ax
            dydt3 = self.gravity      # dvydt = ay

            dydt = np.array([dydt0, dydt1, 
                             dydt2, dydt3])
            return dydt

        #yout = self.euler(dt, yin, derive)
        yout = self.rk2(dt, yin, derive)

        self.posx = yout[0]
        self.posy = yout[1]
        self.velx = yout[2]
        self.vely = yout[3]

        return

    def euler(self, dt, yin, derive):
        """
        Do an Euler update 
        """
        return yout

    def rk2(self, dt, yin, derive):
        """
        Do one RK2 update with dt
        """
        dydt  = derive(yin)
        ystar = yin+dt*dydt
        dydt  = derive(ystar)
        yout  = ystar+dt*dydt
        yout  = 0.5*(yin+yout)
        return yout

        
if __name__=='__main__':


    point = PhysicsBody(initPosx     = 0.0,
                        initPosy     = 0.0,
                        initAngle    = 45.0,
                        initVelocity = 1000.0)

    point.add_gravity()

    for n in np.arange(100):
        point.update(dt=0.01)
        print(point.posx,point.posy)



