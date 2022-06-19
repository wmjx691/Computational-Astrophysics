#binary
from simulation_data import N, stell_info, speedupinfomation
from constants import G, pi, msun, mass, au, yr, sec, dis
from physics import update, update_acc, get_angle, speedup_deter, distance, closest
from output import record
import numpy as np
import copy

from inputfile import ephem, set_starinfo


# initial setup
dt       = 1
t_launch = "1977-09-01 00:00"   #lanuch from earth
t_arrive = "1988-08-15 00:00"   #arrive to   Neptune

# 2dt is prepare to do rk4
position  , velcoity  , time_duration   = ephem(t_launch, t_arrive, 0.5*dt)


# output intervals to record trajectory
step     = 0             #start from 0th step
noutput  = 1000
interval = (len(time_duration)/dt)/noutput


# initialize the stars
mass_star   = np.array([  0.330 *mass,   # 0Mercury
                        4.87  *mass,   # 1Venus
                        5.97  *mass,   # 2Earth
                        0.642 *mass,   # 3Mars
                        1898. *mass,   # 4Jupiter
                        568.  *mass,   # 5Saturn
                        86.8  *mass,   # 6Uranus
                        102.  *mass,   # 7Neptune
                        0.0146*mass,   # 8Pluto
                               msun,   # 9sun
])

radius_star = np.array([  4879,
                        12104,
                        12756,
                        6792,
                        142984,
                        120536,
                        51118,
                        49528,
                        2370,
                        1392700
])
radius_star = radius_star/2

dis_stars   = np.array([57.9,108.2,149.6,227.9,778.6,1433.5,2872.5,4495.1,5906.4,0.])
dis_stars   = dis_stars*dis




# calculate sphere of influence of each planet(soi)
rsoi  =  np.zeros(N)

for i in np.arange(0, N):
    # print(mass_star[0])
    rsoi[i] = dis_stars[i]*(mass_star[i]/mass_star[9])**(2/5)




# stars array has stars[i]:(idindex, mass, x, y)
stars   = set_starinfo(mass_star, position  , velcoity, rsoi)

#output each stars info
first=True
planet=True


# # #create the record .txt of each planet
# position_o, velcoity_o, time_duration_o = ephem(t_launch, t_arrive, dt)
# stars_o = set_starinfo(mass_star, position_o, velcoity_o, rsoi)
# #record each planet in solar system
# record(stars_o, time_duration_o, first, planet)



# 0Mercury
# 1Venus
# 2Earth
# 3Mars
# 4Jupiter
# 5Saturn
# 6Uranus
# 7Neptune
# 8Pluto
# 9sun 

# initialize the stellite
mass_stell    = 721.9    #kg
launch_vel    = 38.621*sec     #km/s->km/day
launch_pos    = 6371. + 1000.   #km


#set info of speed up or not
entersoi_ox    = False
missionindex   = 4       #launch from Earth to Jupiter
speedupornot   = False
speeddirection = np.array([0., 0., 0.])
# speed up unit is 1 day
unit_dt_speed  = 1  #min
dt_ratio       = unit_dt_speed/dt

speedupinfo    = speedupinfomation(entersoi_ox, missionindex, speedupornot, speeddirection, dt_ratio)



# need to get the angle of earth vel vector 
vel_vector_transform  = get_angle(stars[2].vx[0], stars[2].vy[0], stars[2].vz[0])
launch_vel_vector = launch_vel * vel_vector_transform

pos_vector_transform  = get_angle(stars[2].x[0], stars[2].y[0], stars[2].z[0])
launch_pos_vector = launch_pos * pos_vector_transform

stell_setup = np.array([stars[2].x[0] + launch_pos_vector[0],   #x
                        stars[2].y[0] + launch_pos_vector[1],   #y
                        stars[2].z[0] + launch_pos_vector[2],   #z
                        launch_vel_vector[0],  #vx
                        launch_vel_vector[1],  #vy
                        launch_vel_vector[2]   #vz
])



stell = stell_info(         10,  #idindex
                    mass_stell,  #mass
                stell_setup[0],   #x
                stell_setup[1],   #y
                stell_setup[2],   #z
                stell_setup[3],  #vx
                stell_setup[4],  #vy
                stell_setup[5],  #vz
                0,    #ax
                0,    #ay
                0     #az
)

stell_acc = update_acc(step, stell_setup, stell, stars, speedupinfo)

stell.ax  = stell_acc[0]
stell.ay  = stell_acc[1]
stell.az  = stell_acc[2]


# # # test for the earth 2d model
# r        = au
# period   = (((4.0*pi**2)/(G*(msun+stars[2].mass)))*r**3)**0.5
# force    = (G*stars[2].mass*msun)/(r)**2

# stell = stell_info(         10,  #idindex
#                     stars[2].mass,  #mass
#                 r,   #x
#                 0,   #y
#                 0,   #z
#                 0,   #vx
#                 (2.0*pi*r/period),  #vy
#                 0,  #vz
#                 -force/stars[2].mass,    #ax
#                 0,    #ay
#                 0     #az
# )


planet = False

distance0 = distance(stell, stars, speedupinfo.idindex, step)

recordint = 0

#start the simulation
for time in time_duration:
    #prevent the update run over download planet data
    if (step == 7998/dt):
        break

    #rk4 will use the data in the interval of dt, so these step will be skip
    if ((time-2443376.5)%1. != 0):
        step = step+1
        continue

    closest_result, distance0  = closest(stell, stars, speedupinfo, step, distance0)
    speedupinfo                = speedup_deter(stell, stars, speedupinfo, step)

    # check the inital set
    if ((closest_result == True) and (speedupinfo.soi_ox == False)):
        # print("inital launch parameter has problem!!!")
        # print("speed is not enough to go to outer space!!!")
        stell        = update(stell, stars, step, dt, speedupinfo)
        first        = record(stell, time, first, planet)
        step         = step+1

        continue


    # have not arrived target planet
    elif ((closest_result == False) and (speedupinfo.soi_ox == False) ):
        stell        = update(stell, stars, step, dt, speedupinfo)     
        first        = record(stell, time, first, planet)
        step         = step+1

        print("x:", stell.x , "y:", stell.y, "z:", stell.z)

        continue


    # arrive target planet
    elif (speedupinfo.soi_ox == True):
        if (recordint == 0):
            # can back to here
            # step_finddir   = step
            # stell_finddir  = stell

            #the copy problem in happened here
            step_finddir   = copy.copy(step)
            stell_finddir  = copy.copy(stell)
            
            print("find the direction:")
            # first to determine the direction of speed up
            while((closest_result == False)):
                if ((time_duration[step_finddir]-2443376.5)%1. != 0):
                    step_finddir = step_finddir+1
                    continue

                olddis = distance(stell_finddir, stars, speedupinfo.idindex, step_finddir)
                stell_finddir  = update(stell_finddir, stars, step_finddir, dt, speedupinfo)

                step_finddir   = step_finddir+1
                closest_result, olddis  =  closest(stell_finddir, stars, speedupinfo, step_finddir, olddis)

                print("time:",time_duration[step_finddir],"x:", stell_finddir.x , "y:", stell_finddir.y, "z:", stell_finddir.z)

                if (closest_result):
                    step_endsoi = step_finddir

                    x       = stell_finddir.x - stars[speedupinfo.idindex].x[step_endsoi]
                    y       = stell_finddir.y - stars[speedupinfo.idindex].y[step_endsoi]
                    z       = stell_finddir.z - stars[speedupinfo.idindex].z[step_endsoi]

                    speedupinfo.shortest = distance(stell_finddir, stars, speedupinfo.idindex, step_endsoi)

                    speedup_direction     = get_angle(x ,y, z)
                    speedupinfo.direction = speedup_direction
                    
                    break
            print("direction is find!")

            #reset, back to the start point of soi
            speedupinfo.speed_ox = False
            stell  = update(stell, stars, step, dt, speedupinfo)

            print("back to :",time_duration[step],"with no speed")
            print("x:", stell.x , "y:", stell.y, "z:", stell.z)

            closest_result = False

            #compare no speed or speed
            # step_test   = copy.copy(step)
            # stell_test  = copy.copy(stell)

            # speedupinfo.speed_ox = True
            # stell   = update(stell, stars, step, dt, speedupinfo)

            # print("back to :",time_duration[step],"with speed 1 day")
            # print("x:", stell.x , "y:", stell.y, "z:", stell.z)

            first   = record(stell, time, first, planet)
            step    = step+1
            recordint  = recordint+1


        else:
            stell   = update(stell, stars, step, dt, speedupinfo)
            first   = record(stell, time, first, planet)
            step    = step+1


        # # from the soi start point start to jet
        # for step in range(step_startsoi, step_endsoi):

        #     # start to turn on engine 1, 2 ,3,....
        #     for speedtime in range(1, step_endsoi-step_startsoi):

        #         #step return!!!

        #         nowspeed_t = 1
        #         speedupinfo.speed_ox = True
        #         #trajectory of each simulation (with speed up)
        #         while (nowspeed_t <= speedtime):
        #             if ((time_duration[step]-2443376.5)%1. != 0):
        #                 step = step+1
        #                 continue

        #             olddis = distance(stell, stars, speedupinfo.idindex, step)
        #             stell  = update(stell, stars, step, dt, speedupinfo)

        #             nowspeed_t = nowspeed_t+1
        #             step       = step+1

        #             closest_result, olddis  =  closest(stell, stars, speedupinfo, step, olddis)

        #         # record!!
        #         if (closest_result): 
        #             break

        #         speedupinfo.speed_ox = False
        #         #trajectory of each simulation (no engine)
        #         while (closest_result == False):
        #             if ((time_duration[step]-2443376.5)%1. != 0):
        #                 step = step+1
        #                 continue

        #             olddis = distance(stell, stars, speedupinfo.idindex, step)
        #             stell  = update(stell, stars, step, dt, speedupinfo)

        #             step       = step+1

        #             closest_result, olddis  =  closest(stell, stars, speedupinfo, step, olddis)


    # step         = step+1


print("Done for simulation!")