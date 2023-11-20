
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
        sx : Location of source in x (meters)
        sy : Location of source in y (meters)
    output:
        UUm: Acoustic pressure vector (nx,ny) at time n+1
    dependencies:
        None
    written by Jeff Shragge, jshragge@mines.edu, 10/2019
    '''
    ## . . Get dimensions of wavefield
    nx,ny = np.size(UUo,0),np.size(UUo,1)

    ## . . Define Courant numbers (squared)
    dtdx2,dtdy2 = (dt/dx)**2,(dt/dy)**2

    ## Source location
    isx,isy = int(sx/dx),int(sy/dy) ## . . Force to be integer

    ## . . Inject wavelet
    UUo[isx,isy] += dt*dt*F[it]

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