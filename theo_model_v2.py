import numpy as np
import matplotlib.pyplot as plt

## . . Define spatial grid
Lx,Ly = 2000,1200           # . . lengths (m)
nx,ny = 401,241             # . . Number of points in discretization
dx,dy = Lx/(nx-1),Ly/(ny-1) # . . Discretization intervals

## . . Create velocity function
v = np.zeros((nx,ny))+1500
v[:,50:]=2000


import numpy as np
def awe_2d_explicit_solver_heterogeneous_8th_order(UUo,UUm,dx,dy,dt,v,F,it,sx,sy):
    '''Set up eight-order solver of the acoustic wave equation
    usage: U=awe_2d_explicit_solver_heterogeneous_8th_order(UUo,UUm,dx,dy,dt,v,F,it,sx,sy):
    input:
        UUo: Acoustic pressure vector (nx,ny) at time step n
        UUm: Acoustic pressure vector (nx,ny) at time step n-1
        dx : Spatial sampling in x
        dy : Spatial sampling in y
        dt : Temporal sampling
        v  : Heterogeneous propagation velocity (nx,ny)
        F  : Forcing function (nt)
        it : Time index
        sx : Location of source in x (meters). In this case, sx should be an array.
        sy : Location of source in y (meters). In this case, sx should be an array.
    output:
        UUm: Acoustic pressure vector (nx,ny) at time n+1
    dependencies:
        None
    written by Jeff Shragge, jshragge@mines.edu, 10/2019
    extended to multiple source case by Shenyao Jin, shenyaojin@mines.edu, 11/2023
    '''
    ## . . Get dimensions of wavefield
    nx,ny = np.size(UUo,0),np.size(UUo,1)

    ## . . Define Courant numbers (squared)
    dtdx2,dtdy2 = (dt/dx)**2,(dt/dy)**2

    ## Source location (iteration)

    if isinstance(sx, (int, float)):
        isx,isy = int(sx/dx),int(sy/dy)
        # inject wavelet
        UUo[isx,isy] += dt*dt*F[it]
    else:
        for iter in range(len(sx)):
            isx = int(sx[iter]/dx)
            isy = int(sy[iter]/dy)
            UUo[isx, isy] += dt * dt * F[it]
    # isx,isy = int(sx/dx),int(sy/dy) ## . . Force to be integer


    ## . . Update solution
    UUm[4:nx-4,4:ny-4] =  2*UUo[4:nx-4,4:ny-4]-UUm[4:nx-4,4:ny-4]+dtdx2*v[4:nx-4,4:ny-4]**2*(
            -1/560 *UUo[0:nx-8,4:ny-4]
            +8/315 *UUo[1:nx-7,4:ny-4]
            -1/5   *UUo[2:nx-6,4:ny-4]
            +8/5   *UUo[3:nx-5,4:ny-4]
            -205/72*UUo[4:nx-4,4:ny-4]
            +8/5   *UUo[5:nx-3,4:ny-4]
            -1/5   *UUo[6:nx-2,4:ny-4]
            +8/315 *UUo[7:nx-1,4:ny-4]
            -1/560 *UUo[8:nx  ,4:ny-4])+dtdy2*v[4:nx-4,4:ny-4]**2*(
                                  -1/560 *UUo[4:nx-4,0:ny-8]
                                  +8/315 *UUo[4:nx-4,1:ny-7]
                                  -1/5   *UUo[4:nx-4,2:ny-6]
                                  +8/5   *UUo[4:nx-4,3:ny-5]
                                  -205/72*UUo[4:nx-4,4:ny-4]
                                  +8/5   *UUo[4:nx-4,5:ny-3]
                                  -1/5   *UUo[4:nx-4,6:ny-2]
                                  +8/315 *UUo[4:nx-4,7:ny-1]
                                  -1/560 *UUo[4:nx-4,8:ny  ])

    return UUm ## . . Return updated wavefield at time step n+1

# also, include the reverse pde solver
def awe_2d_heterogeneous_8th_order_data_time_reverse(UUo,UUm,dx,dy,dt,v,D,it,ry):
    '''Set up eight-order solver of the acoustic wave equation
    usage: U=awe_2d_heterogeneous_8th_order_data_time_reverse(UUo,UUm,dx,dy,dt,v,D,it,ry):
    input:
        Uo: Acoustic pressure vector (nx,ny) at time step n
        Um: Acoustic pressure vector (nx,ny) at time step n-1
        dx: Spatial sampling in x
        dy: Spatial sampling in y
        dt: Temporal sampling
        v : Heterogeneous propagation velocity (nx,ny)
        D : Data to be time reversed
        it: Time index
        ry: Injection depth (meters)
    output:
        Um: Acoustic pressure vector (nx,ny) at time n+1
    dependencies:
        None
    written by Jeff Shragge, jshragge@mines.edu, 10/2019
    '''
    ## . . Get dimensions of wavefield
    nx,ny = np.size(UUo,0),np.size(UUo,1)

    ## . . Define Courant numbers (squared)
    dtdx2 = (dt/dx)**2
    dtdy2 = (dt/dy)**2

    ## Source location
    iry = int(ry/dy) ## . . Force to be integer

    ## . . Inject wavelet
    UUo[:,iry] += dt*dt*D[:,it]

    ## . . Update solution
    UUm[4:nx-4,4:ny-4] =  2*UUo[4:nx-4,4:ny-4]-UUm[4:nx-4,4:ny-4]+dtdx2*v[4:nx-4,4:ny-4]**2*(
            -1/560 *UUo[0:nx-8,4:ny-4]
            +8/315 *UUo[1:nx-7,4:ny-4]
            -1/5   *UUo[2:nx-6,4:ny-4]
            +8/5   *UUo[3:nx-5,4:ny-4]
            -205/72*UUo[4:nx-4,4:ny-4]
            +8/5   *UUo[5:nx-3,4:ny-4]
            -1/5   *UUo[6:nx-2,4:ny-4]
            +8/315 *UUo[7:nx-1,4:ny-4]
            -1/560 *UUo[8:nx  ,4:ny-4])+dtdy2*v[4:nx-4,4:ny-4]**2*(
                                  -1/560 *UUo[4:nx-4,0:ny-8]
                                  +8/315 *UUo[4:nx-4,1:ny-7]
                                  -1/5   *UUo[4:nx-4,2:ny-6]
                                  +8/5   *UUo[4:nx-4,3:ny-5]
                                  -205/72*UUo[4:nx-4,4:ny-4]
                                  +8/5   *UUo[4:nx-4,5:ny-3]
                                  -1/5   *UUo[4:nx-4,6:ny-2]
                                  +8/315 *UUo[4:nx-4,7:ny-1]
                                  -1/560 *UUo[4:nx-4,8:ny  ])

    return UUm ## . . Return updated wavefield at time step n+1

## . . Init wavefields on spatial grid
UUo = np.zeros((nx,ny))
UUm = np.zeros((nx,ny))

## Time stepping parameters
CC = 0.5                     # . . Courant #
nt = 800                     # . . Number of time steps
dt = CC*dx/np.max(v)         # . . Define dt based on Courant
t  = np.linspace(0,nt*dt,nt) # . . Time lin
t0 = 0.05                    # . . Wavelet shift

## . . Define forcing function
ss=0.01                      # . . sigma for Ricker wavelet
F = (1-((t-t0)/ss)**2)*np.exp(-(t-t0)**2/(2*ss**2))
shot_num = 10
sx = np.linspace(0, nx * dx, shot_num + 2)[1:-1]
sy = 25

img = np.zeros((nx, ny))
for iter in range(shot_num):
    fff = np.zeros((nx,ny,nt))
    for it in range(nt):
        tmp = awe_2d_explicit_solver_heterogeneous_8th_order(UUo,UUm,dx,dy,dt,v,F,it,sx[iter],sy) #calc solution at n+1
        fff[:,:,it]=tmp  ## save solution vector
        UUm=UUo          ## move solution at n to n-1 to prepare for next iteration
        UUo=tmp          ## move solution at n+1 to n to prepare for next iteration

    ## . . Create velocity function
    v = np.zeros((nx,ny))+1500
    v[:,50:]=2000

    ## . . Init wavefields on spatial grid
    UUo = np.zeros((nx,ny))
    UUm = np.zeros((nx,ny))

    ## Time stepping parameters
    CC = 0.5                     # . . Courant #
    nt = 800                     # . . Number of time steps
    dt = CC*dx/np.max(v)         # . . Define dt based on Courant
    t  = np.linspace(0,nt*dt,nt) # . . Time lin
    t0 = 0.05                    # . . Wavelet shift

    ## . . Define forcing function
    ss=0.01                      # . . sigma for Ricker wavelet

    ## . . Get the data from fff above
    data = fff[:,5,:]

    ## . . Define source location
    ry=25                        # . . receiver injection location in y (in physical m units)

    ## . . Total Solution space
    mmm = np.zeros((nx,ny,nt))

    ## . . Iterate over solution
    ## . . Note time reversal
    for it in range(nt-1,0,-1):
        tmp = awe_2d_heterogeneous_8th_order_data_time_reverse(UUo,UUm,dx,dy,dt,v,data,it,ry) #calc solution at n+1
        mmm[:,:,nt-1-it]=tmp  ## save solution vector
        UUm=UUo          ## move solution at n to n-1 to prepare for next iteration
        UUo=tmp          ## move solution at n+1 to n to prepare for next iteration

    img += np.sum(fff*mmm[:,:,::-1], axis=2)

cx = np.array([-1,1])
plt.imshow(img.T, aspect='auto',cmap='jet',extent=[0,Lx,Ly,0])
plt.clim(cx * 5e-16)
plt.colorbar()
plt.show()