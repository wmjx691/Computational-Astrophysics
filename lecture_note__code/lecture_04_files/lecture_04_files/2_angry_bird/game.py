import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.widgets import Slider, Button
from PIL import Image

from character import AngryBird, HappyPig
from physics import PhysicsBody

class AngryBirdGame:

    def __init__(self, angle, velocity, gravity, dt,
                 xmax=1200, ymax=600, buf=90):

        # set default vales
        self.angle    = angle
        self.velocity = velocity
        self.gravity  = gravity
        self.dt       = dt
        self.xmax     = xmax
        self.ymax     = ymax
        self.run      = False

        # Bird
        size = 90
        self.size  = size
        self.hsize = int(size/2)
        bird = AngryBird(name = "bird",
                         img = "./images/bird.png")
        bird.setup_image(size=size)

        # bird has a physics body
        body = PhysicsBody(initPosx     = 0.0,
                           initPosy     = 0.0,
                           initAngle    = angle,
                           initVelocity = velocity)
        body.add_gravity(gravity=gravity)
        bird.body = body
        self.bird = bird

        # Pig
        pig = HappyPig(name = "pig",
                       img = "./images/pig.jpg")
        pig.setup_image(size=size)

        body = PhysicsBody(initPosx     = xmax -100,
                           initPosy     = 0.0,
                           initAngle    = 0.0,
                           initVelocity = 0.0)
        body.add_gravity(gravity=0.0)
        pig.body = body
        self.pig = pig

        self.setup_frame(xmax,ymax,buf)
        
        return


    def setup_frame(self, xmax, ymax, buf=50):

        fig, ax = plt.subplots()
        plt.subplots_adjust(left=0.1, bottom=0.25)

        # draw the initial predicited tgrajectory
        tt,xx,yy = self.get_predicted_trajectory(self.angle, self.velocity)
        traj, = ax.plot(xx[20:-200],yy[20:-200],color='0.8',animated=False)
        self.traj = traj

        # set axis
        ax.set_xlim([-buf, xmax+buf])
        ax.set_ylim([-buf, ymax+buf])
        ax.axis('off')

        # draw bird/pig     #assume size =90
        size  = self.size
        hsize = self.hsize

        im_bird = ax.imshow(self.bird.image,extent = (-hsize,hsize,-hsize,hsize))
        im_pig  = ax.imshow(self.pig.image,extent  = (xmax-100,xmax-100+size,-hsize,hsize))
        self.bird.imshow = im_bird
        self.pig.imshow  = im_pig

        ### setup widgets

        # sliders
        axcolor = 'lightgoldenrodyellow'
        ax_ang  = plt.axes([0.2, 0.1,  0.65, 0.03], facecolor=axcolor)
        ax_vel  = plt.axes([0.2, 0.15, 0.65, 0.03], facecolor=axcolor)

        slider_ang = Slider(ax_ang, 'Angle', 0.0, 90.0, valinit=self.angle)
        slider_vel = Slider(ax_vel, 'Velocity', 300.0, 2000.0, valinit=self.velocity, valstep=20.0)

        self.slider_ang = slider_ang
        self.slider_vel = slider_vel

        slider_ang.on_changed(self.update_trajectory)
        slider_vel.on_changed(self.update_trajectory)

        # buttons
        resetax = plt.axes([0.75, 0.025, 0.1, 0.04])
        button_reset = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

        shootax = plt.axes([0.2,0.025, 0.1,0.04])
        button_shoot = Button(shootax, 'Shoot', color=axcolor, hovercolor='0.975')

        self.button_reset = button_reset
        self.button_shoot = button_shoot

        def reset(event):
            self.run = False
            self.slider_ang.reset()
            self.slider_vel.reset()
            self.init_bird()
            return

        button_reset.on_clicked(reset)

        def shoot(event):
            self.run = True
            bird = self.bird
            bird.body.reset(
                    initPosx = 0.0,
                    initPosy = 0.0,
                    initAngle = self.slider_ang.val,
                    initVelocity = self.slider_vel.val
                    )
            return

        button_shoot.on_clicked(shoot)

        self.ax  = ax
        self.fig = fig
        return

    def get_predicted_trajectory(self, angle, velocity):
        angle = angle / 180.0 * np.pi
        velx  = velocity * np.cos(angle)
        vely  = velocity * np.sin(angle)

        xx = np.arange(0,self.xmax,2.0)
        tt = xx / velx     
        yy = vely*tt + 0.5*self.gravity*tt**2

        return tt,xx,yy

    def update_trajectory(self, val):
        angle    = self.slider_ang.val
        velocity = self.slider_vel.val

        if not self.run:
            tt,xx,yy = self.get_predicted_trajectory(angle, velocity)
            self.traj.set_xdata(xx[20:-200])
            self.traj.set_ydata(yy[20:-200])
        return

    def update(self, i):

        dt_out   = 0.1     # sec
        dt       = self.dt # sec
        nstep    = max(int(dt_out/dt),1)
        # redraw bird's location
        if self.run:

            self.bird.body.update(dt=self.dt)

            if (i%nstep == 0):
                self.draw(self.bird)

            self.check_result()

        return

    def draw(self, char):

        hs = self.hsize
        ax = self.ax
        # assume there is a body
        body = char.body
        empty   = np.add.outer(range(4),range(4))*np.nan
        im_char = char.imshow
        im_char.set_data(empty)
        im_char = ax.imshow(char.image, extent = 
                                (-hs+body.posx,hs+body.posx,-hs+body.posy,hs+body.posy))
        char.imshow = im_char
        return

    def init_bird(self):
        # should have a better way
        ax = self.ax
        im_bird = self.bird.imshow
        empty = np.add.outer(range(4),range(4))*np.nan
        im_bird.set_data(empty)
        im_bird = ax.imshow(self.bird.image,extent = (-45,45,-45,45))
        self.bird.imshow = im_bird
        return

    def check_result(self):

        bird = self.bird
        pig  = self.pig
        hs   = self.hsize

        if (bird.body.posy < 0.0):
            # ended
            posx = bird.body.posx
            win_zone = [pig.body.posx-hs,
                        pig.body.posx+hs]

            if ((posx >= win_zone[0]) and
                (posx <= win_zone[1])):
                print("You Win !!!")
            else:
                print("Game Over!")
            self.run = False

        return

    def play(self):

        fig = self.fig

        animate = animation.FuncAnimation(fig, self.update,
                                          frames=100, interval=1, 
                                          init_func=self.init_bird,
                                          blit=False)

        plt.show()
        return

if __name__=='__main__':

    default_angle    = 60.0     # degree
    default_velocity = 1000.0   # cm/s
    default_gravity  = -980.0   # cm/s/s
    default_dt       = 0.1      # sec

 
    game = AngryBirdGame(angle    = default_angle,
                         velocity = default_velocity,
                         gravity  = default_gravity,
                         dt       = default_dt,
                         xmax     = 1200,
                         ymax     = 600)

    game.play()


