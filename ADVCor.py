import numpy as np
from scipy import ndimage
import scipy.integrate
from relax import relax2d,relax3d

def GalChen2D(field1, field2, delta_t,dx,missing=999.):
    
    """This function finds a crude two-dimensional Gal-Chen first guess for the
       horizontal pattern translation components.
    """
    
    temp_field1 = np.copy(field1)
    temp_field2 = np.copy(field2)
    
    temp_field1[temp_field1 == missing] = np.nan
    temp_field2[temp_field2 == missing] = np.nan
    
    rdx = 1./dx
    
    dRdt = (temp_field2[1:field2.shape[0]-1,1:field2.shape[1]-1]-temp_field1[1:field2.shape[0]-1,1:field2.shape[1]-1])/delta_t
    
    dRdy = 0.25*rdx*(temp_field1[2:field1.shape[0],1:field1.shape[1]-1]-temp_field1[0:field1.shape[0]-2,1:field1.shape[1]-1]+
                     temp_field2[2:field2.shape[0],1:field2.shape[1]-1]-temp_field2[0:field2.shape[0]-2,1:field2.shape[1]-1])
    
    dRdx = 0.25*rdx*(temp_field1[1:field1.shape[0]-1,2:field1.shape[1]]-temp_field1[1:field1.shape[0]-1,0:field1.shape[1]-2]+
                     temp_field2[1:field2.shape[0]-1,2:field2.shape[1]]-temp_field2[1:field2.shape[0]-1,0:field2.shape[1]-2])
    
    
    a = dx * dx * dRdt * dRdx
    b = dx * dx * dRdx * dRdx
    c = dx * dx * dRdx * dRdy
    d = dx * dx * dRdt * dRdy
    e = dx * dx * dRdy * dRdy
    
    suma = np.nansum(a)
    sumb = np.nansum(b)
    sumc = np.nansum(c)
    sumd = np.nansum(d)
    sume = np.nansum(e)
    
    temp_u = (suma*sume - sumc*sumd)/(sumc*sumc - sumb*sume)
    temp_v = (sumb*sumd - sumc*suma)/(sumc*sumc - sumb*sume)
    
    u = np.ones(field1.shape)*temp_u
    v = np.ones(field1.shape)*temp_v

    return u,v


def GalChen3D(field1, field2, delta_t, dx, missing=999.):
    
    """This function finds a crude three-dimensional Gal-Chen first guess 
       for the horizontal and vertical pattern translation components.
    """
    
    temp_field1 = np.copy(field1)
    temp_field2 = np.copy(field2)
    
    temp_field1[temp_field1 == missing] = np.nan
    temp_field2[temp_field2 == missing] = np.nan
    
    rdx = 1./dx
    
    dRdt = (temp_field2[1:field2.shape[0]-1,1:field2.shape[1]-1,1:field2.shape[2]-1]-
            temp_field1[1:field1.shape[0]-1,1:field1.shape[1]-1,1:field1.shape[2]-1])/delta_t
    
    dRdz = 0.25*rdx*(temp_field1[2:,1:field1.shape[1]-1,1:field1.shape[2]-1]-temp_field1[0:field1.shape[0]-2,1:field1.shape[1]-1,1:field1.shape[2]-1]+
                     temp_field2[2:,1:field2.shape[1]-1,1:field2.shape[2]-1]-temp_field2[0:field2.shape[0]-2,1:field2.shape[1]-1,1:field2.shape[2]-1])
    
    dRdy = 0.25*rdx*(temp_field1[1:field1.shape[0]-1,2:,1:field1.shape[2]-1]-temp_field1[1:field1.shape[0]-1,0:field1.shape[1]-2,1:field1.shape[2]-1]+
                     temp_field2[1:field2.shape[0]-1,2:,1:field2.shape[2]-1]-temp_field2[1:field2.shape[0]-1,0:field2.shape[1]-2,1:field2.shape[2]-1])
    
    dRdx = 0.25*rdx*(temp_field1[1:field1.shape[0]-1,1:field1.shape[1]-1,2:]-temp_field1[1:field1.shape[0]-1,1:field1.shape[2]-1,0:field1.shape[2]-2]+
                     temp_field2[1:field2.shape[0]-1,1:field2.shape[1]-1,2:]-temp_field2[1:field2.shape[0]-1,1:field2.shape[2]-1,0:field2.shape[2]-2])
    
    a = dx * dx * dRdt * dRdx
    b = dx * dx * dRdx * dRdx
    c = dx * dx * dRdx * dRdy
    d = dx * dx * dRdx * dRdz
    e = dx * dx * dRdt * dRdy
    f = dx * dx * dRdy * dRdy
    g = dx * dx * dRdy * dRdz
    h = dx * dx * dRdt * dRdz
    l = dx * dx * dRdz * dRdz
    
    suma = np.nansum(a)
    sumb = np.nansum(b)
    sumc = np.nansum(c)
    sumd = np.nansum(d)
    sume = np.nansum(e)
    sumf = np.nansum(f)
    sumg = np.nansum(g)
    sumh = np.nansum(h)
    suml = np.nansum(l)
    
    temp_u = ((sumc*(sumg*sumh-suml*sume)+
              sumd*(sumg*sume-sumf*sumh)+
              suma*(sumf*suml-sumg*sumg))/
              (sumb*(sumg*sumg-sumf*suml)+
               sumc*sumc*suml-
               2.*sumc*sumd*sumg+
               sumd*sumd*sumf))
    
    temp_v = -((sumb*(sumg*sumh-suml*sume)+
               sumd*sumd*sume+
               suma*sumc*suml-
               sumd*(sumc*sumh+suma*sumg))/
              (sumb*(sumg*sumg-sumf*suml)+
               sumc*sumc*suml-
               2.*sumc*sumd*sumg+
               sumd*sumd*sumf))
    
    temp_w = ((sumb*(sumf*sumh-sumg*sume)+
              sumd*(sumc*sume-suma*sumf)-
              sumc*sumc*sumh+
              suma*sumc*sumg)/
              (sumb*(sumg*sumg-sumf*suml)+
               sumc*sumc*suml-
               2.*sumc*sumd*sumg+
               sumd*sumd*sumf))
    
    
    u = np.ones(field1.shape)*temp_u
    v = np.ones(field1.shape)*temp_v
    w = np.ones(field1.shape)*temp_w

    return u,v,w    


def ADV2D(field1, field2, first_U, first_V, dx, dy, bigT, nt, dt, beta,relax=1,under=1,itermax=20000,
          itermainmax=100,tol=0.001,tol2=0.01,missing=999.,verbose=True):
    
    
    """This function performs the 2-D advection correction as described in Shapiro
       et al. 2010.
    """

    u = np.copy(first_U)
    v = np.copy(first_V)
    
    rdx = 1./dx
    rdy = 1./dy
    rdt = 1./dt
    
    x = np.arange(field1.shape[0])*dx
    y = np.arange(field2.shape[1])*dy
    
    # Check to make sure the timing works out for bigT, nt, and dt
    
    if int(bigT/dt) != (nt-1):
        raise ValueError('dt*(nt-1) does not equal bigT')
    
    # Big loops for the main iteration. This is the biggest inefficiency with the
    # Python version of this procedure
    if verbose:
        print('Entering main iterative loop')
    
    div = np.ones(field1.shape)*np.nan
    
    div[1:div.shape[0]-1,1:div.shape[1]-1] = (0.5*rdx*(u[1:u.shape[0]-1,2:u.shape[1]] - u[1:u.shape[0]-1,0:u.shape[1]-2]) +
                                              0.5*rdy*(v[2:v.shape[0],1:v.shape[1]-1] - v[0:v.shape[0]-2,1:v.shape[1]-1]))
    div[:,0] = div[:,1]
    div[:,-1] = div[:,-2]
    
    div[0,:] = div[1,:]
    div[-1,:] = div[-2,:]
    
    
    for itermain in range(1,itermainmax+1):
        
        if verbose:
            print('itermain = ' + str(itermain))
        
        forwardref = np.ones((field1.shape[0],field1.shape[1],nt))*np.nan
        backwardref = np.ones((field1.shape[0],field1.shape[1],nt))*np.nan
        divtraj = np.ones((field1.shape[0],field1.shape[1],nt,nt))*np.nan
        
        
        forwardref[:,:,0] = field1
        forwardref[:,:,-1] = field2
        backwardref[:,:,0] = field1
        backwardref[:,:,-1] = field2
        
        divtraj[:,:,0,0] = div[:,:]
        divtraj[:,:,-1,-1] = div[:,:]
        
        xtemp_f = (np.ones(field1.shape)*x[None,:])
        ytemp_f = (np.ones(field1.shape)*y[:,None])
        
        xtemp_b = (np.ones(field1.shape)*x[None,:])
        ytemp_b = (np.ones(field1.shape)*y[:,None])
        
        for traj in range(1,nt-1):
            
            divtraj[:,:,traj,traj] = div[:,:]
            # First calculate forward trajectories
            
            varout = ndimage.map_coordinates(u,[ytemp_f.ravel()*rdy,xtemp_f.ravel()*rdx],order=1,cval=np.nan)
            k1x = dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(v,[ytemp_f.ravel()*rdy,xtemp_f.ravel()*rdx],order=1,cval=np.nan)
            k1y = dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(u,[(ytemp_f+(0.5*k1y)).ravel()*rdy,(xtemp_f+(0.5*k1x)).ravel()*rdx],order=1,cval=np.nan)
            k2x = dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(v,[(ytemp_f+(0.5*k1y)).ravel()*rdy,(xtemp_f+(0.5*k1x)).ravel()*rdx],order=1,cval=np.nan)
            k2y = dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(u,[(ytemp_f+(0.5*k2y)).ravel()*rdy,(xtemp_f+(0.5*k2x)).ravel()*rdx],order=1,cval=np.nan)
            k3x = dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(v,[(ytemp_f+(0.5*k2y)).ravel()*rdy,(xtemp_f+(0.5*k2x)).ravel()*rdx],order=1,cval=np.nan)
            k3y = dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(u,[(ytemp_f+k3y).ravel()*rdy,(xtemp_f+k3x).ravel()*rdx],order=1,cval=np.nan)
            k4x = dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(v,[(ytemp_f+k3y).ravel()*rdy,(xtemp_f+k3x).ravel()*rdx],order=1,cval=np.nan)
            k4y = dt * np.copy(varout.reshape(field1.shape))
            
            xtemp_f = xtemp_f + (k1x + 2.*k2x + 2.*k3x + k4x)/6.
            ytemp_f = ytemp_f + (k1y + 2.*k2y + 2.*k3y + k4y)/6.
            
            ind = np.arange(1,nt-traj)
            varout = ndimage.map_coordinates(div,[ytemp_f.ravel()*rdy,xtemp_f.ravel()*rdx],order=1,cval=np.nan)
            divtraj[:,:,ind,ind+traj] = np.array([(np.copy(varout.reshape(field1.shape))).T]*(nt-1-traj)).T
            
            varout = ndimage.map_coordinates(field2,[ytemp_f.ravel()*rdy,xtemp_f.ravel()*rdx],order=1,cval=np.nan)
            forwardref[:,:,nt-1-traj] = np.copy(varout.reshape(field1.shape))
            
            # Now do backward trajectories
            
            varout = ndimage.map_coordinates(u,[ytemp_b.ravel()*rdy,xtemp_b.ravel()*rdx],order=1,cval=np.nan)
            k1x = -dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(v,[ytemp_b.ravel()*rdy,xtemp_b.ravel()*rdx],order=1,cval=np.nan)
            k1y = -dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(u,[(ytemp_b+(0.5*k1y)).ravel()*rdy,(xtemp_b+(0.5*k1x)).ravel()*rdx],order=1,cval=np.nan)
            k2x = -dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(v,[(ytemp_b+(0.5*k1y)).ravel()*rdy,(xtemp_b+(0.5*k1x)).ravel()*rdx],order=1,cval=np.nan)
            k2y = -dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(u,[(ytemp_b+(0.5*k2y)).ravel()*rdy,(xtemp_b+(0.5*k2x)).ravel()*rdx],order=1,cval=np.nan)
            k3x = -dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(v,[(ytemp_b+(0.5*k2y)).ravel()*rdy,(xtemp_b+(0.5*k2x)).ravel()*rdx],order=1,cval=np.nan)
            k3y = -dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(u,[(ytemp_b+k3y).ravel()*rdy,(xtemp_b+k3x).ravel()*rdx],order=1,cval=np.nan)
            k4x = -dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(v,[(ytemp_b+k3y).ravel()*rdy,(xtemp_b+k3x).ravel()*rdx],order=1,cval=np.nan)
            k4y = -dt * np.copy(varout.reshape(field1.shape))
            
            xtemp_b = xtemp_b + (k1x + 2.*k2x + 2.*k3x + k4x)/6.
            ytemp_b = ytemp_b + (k1y + 2.*k2y + 2.*k3y + k4y)/6.
            
            ind = np.arange(traj,nt-1)
            varout = ndimage.map_coordinates(div,[ytemp_b.ravel()*rdy,xtemp_b.ravel()*rdx],order=1,cval=np.nan)
            divtraj[:,:,ind,ind-traj] = np.array([(np.copy(varout.reshape(field1.shape))).T]*(nt-1-traj)).T
            
            varout = ndimage.map_coordinates(field1,[ytemp_b.ravel()*rdy,xtemp_b.ravel()*rdx],order=1,cval=np.nan)
            backwardref[:,:,traj] = np.copy(varout.reshape(field1.shape))
        
        # Now you have to calculate the integral of the divergence of u and v
        # along the trajectories
        
        gral1 = np.zeros(divtraj.shape)
        
        gral1[:,:,:,1:] = scipy.integrate.cumtrapz(divtraj,dx=dt,axis=3)

        # and another intergral which is the integral of the exponential of the 
        # previous integral

        gral2 = np.zeros(divtraj.shape)

        gral2[:,:,:,1:] = scipy.integrate.cumtrapz(np.exp(-gral1),dx=dt,axis=3)
        
        # Here we construct the reflectivity field using the trajectories
        # and the integrals we just calculated
        
        ratio = gral2[:,:,np.arange(nt),np.arange(nt)]/gral2[:,:,np.arange(nt),-1]
        
        ref = backwardref + (ratio * (forwardref - backwardref))
        
        ref[:,:,0] = backwardref[:,:,0]
        ref[:,:,-1] = backwardref[:,:,-1]
        
        # Here we caluclate the coefficients for the pde for u and v. This is 
        # similar to the Gal-Chen (1982) technique
        
        dRdt = 0.5 * rdt * (ref[1:ref.shape[0]-1,1:ref.shape[1]-1,2:] - 
                            ref[1:ref.shape[0]-1,1:ref.shape[1]-1,:ref.shape[2]-2])
        
        dRdy = 0.5 * rdy * (ref[2:,1:ref.shape[1]-1,1:ref.shape[2]-1] -
                            ref[:ref.shape[0]-2,1:ref.shape[1]-1,1:ref.shape[2]-1])
        
        dRdx = 0.5 * rdx * (ref[1:ref.shape[0]-1,2:,1:ref.shape[2]-1] -
                            ref[1:ref.shape[0]-1,:ref.shape[1]-2,1:ref.shape[2]-1])
        
        foo = np.where((np.isnan(dRdt)) | (np.isnan(dRdy)) | (np.isnan(dRdx)))
        
        dRdt[foo] = np.nan
        dRdx[foo] = np.nan
        dRdy[foo] = np.nan
        
        a = np.zeros(u.shape)
        b = np.zeros(u.shape)
        c = np.zeros(u.shape)
        d = np.zeros(u.shape)
        e = np.zeros(u.shape)
         
        a[1:ref.shape[0]-1,1:ref.shape[1]-1] = np.nansum(dx * dx * dRdt * dRdx / beta, axis = 2)/(nt-2)
        b[1:ref.shape[0]-1,1:ref.shape[1]-1] = np.nansum(dx * dx * dRdx * dRdx / beta, axis = 2)/(nt-2)
        c[1:ref.shape[0]-1,1:ref.shape[1]-1] = np.nansum(dx * dx * dRdx * dRdy / beta, axis = 2)/(nt-2)
        d[1:ref.shape[0]-1,1:ref.shape[1]-1] = np.nansum(dx * dx * dRdt * dRdy / beta, axis = 2)/(nt-2)
        e[1:ref.shape[0]-1,1:ref.shape[1]-1] = np.nansum(dx * dx * dRdy * dRdy / beta, axis = 2)/(nt-2)
        
        # Now we call the function that performs the successive overrelaxation
        # to solve the PDEs for U and V. I am unaware of any Python packages
        # that have PDE solvers that would work here so this function is a 
        # f2py'ed fortran function
        
        uold = np.copy(u)
        vold = np.copy(v)
        
        u_new, v_new = relax2d(a,b,c,d,e,u,v,itermax,tol,relax)
        
        # Check to see if u, v have converged in the big iteration
        
        dif = np.sqrt((u_new-uold)**2 + (v_new-vold)**2)
        if np.max(dif) < tol2:
            if verbose:
                print("*** Convergence in big iteration loop ***")
            break
        else:
            if verbose:
                print(" ")
                print("itermain = " + str(itermain))
                print("U relaxation solution in center of domain: " + str(u_new[int(u.shape[0]/2),int(u.shape[1]/2)]))
                print("V relaxation solution in center of domain: " + str(v_new[int(u.shape[0]/2),int(u.shape[1]/2)]))
                print("Maximum difference: " + str(np.max(dif)))
        
            u = (1.0 - under) * uold + under * u_new
            v = (1.0 - under) * vold + under * v_new
            
            div[1:div.shape[0]-1,1:div.shape[1]-1] = (0.5*rdx*(u[1:u.shape[0]-1,2:u.shape[1]] - u[1:u.shape[0]-1,0:u.shape[1]-2]) +
                                              0.5*rdy*(v[2:v.shape[0],1:v.shape[1]-1] - v[0:v.shape[0]-2,1:v.shape[1]-1]))
            div[:,0] = div[:,1]
            div[:,-1] = div[:,-2]
    
            div[0,:] = div[1,:]
            div[-1,:] = div[-2,:]
        
    return u, v, ref


def ADV3D(field1, field2, first_U, first_V, first_W, dx, dy, dz, bigT, nt, dt, beta, gamma, eta, nu,
          relax=1,under=1,itermax=20000,itermainmax=100,tol=0.001,tol2=0.01,missing=999.,verbose=True):
    
    """This function performs the 3-D advection correction as described in
       Gebauer et al. 2021 (In prep).
    """
    
    u = np.copy(first_U)
    v = np.copy(first_V)
    w = np.copy(first_W)
    
    rdx = 1./dx
    rdy = 1./dy
    rdz = 1./dz
    rdt = 1./dt
    
    x = np.arange(field1.shape[2])*dx
    y = np.arange(field1.shape[1])*dy
    z = np.arange(field1.shape[0])*dz
    
    # Check to make sure the timing works out for bigT, nt, and dt
    
    if int(bigT/dt) != (nt-1):
        raise ValueError('dt*(nt-1) does not equal bigT')
    
    # Big loops for the main iteration. This is the biggest inefficiency with the
    # Python version of this procedure
    if verbose:
        print('Entering main iterative loop')
    
    div = np.ones(field1.shape)*np.nan
    
    div[1:div.shape[0]-1,1:div.shape[1]-1,1:div.shape[2]-1] = (
        0.5*rdz*(w[2:,1:w.shape[1]-1,1:w.shape[2]-1] - w[:u.shape[0]-2,1:w.shape[1]-1,1:w.shape[2]-1]) +
        0.5*rdy*(v[1:v.shape[0]-1,2:,1:v.shape[2]-1] - v[1:v.shape[0]-1,:v.shape[1]-2,1:v.shape[1]-1]) +
        0.5*rdx*(u[1:u.shape[0]-1,1:u.shape[1]-1,2:] - u[1:v.shape[0]-1,1:v.shape[1]-1,:v.shape[2]-2])
        )
    
    div[:,0,:] = div[:,1,:]
    div[:,-1,:] = div[:,-2,:]
    
    div[:,:,0] = div[:,:,1]
    div[:,:,-1] = div[:,:,-2]
    
    div[0,:,:] = div[1,:,:]
    div[-1,:,:] = div[-2,:,:]
    
    for itermain in range(1,itermainmax+1):
        
        if verbose:
            print('itermain = ' + str(itermain))
        
        forwardref = np.ones((field1.shape[0],field1.shape[1],field1.shape[2],nt))*np.nan
        backwardref = np.ones((field1.shape[0],field1.shape[1],field1.shape[2],nt))*np.nan
        divtraj = np.ones((field1.shape[0],field1.shape[1],field1.shape[2],nt,nt))*np.nan
        
        forwardref[:,:,:,0] = field1
        forwardref[:,:,:,-1] = field2
        backwardref[:,:,:,0] = field1
        backwardref[:,:,:,-1] = field2
        
        divtraj[:,:,:,0,0] = div[:,:,:]
        divtraj[:,:,:,-1,-1] = div[:,:,:]
        
        xtemp_f = (np.ones(field1.shape)*x[None,None,:])
        ytemp_f = (np.ones(field1.shape)*y[None,:,None])
        ztemp_f = (np.ones(field1.shape)*z[:,None,None])
        
        xtemp_b = (np.ones(field1.shape)*x[None,None,:])
        ytemp_b = (np.ones(field1.shape)*y[None,:,None])
        ztemp_b = (np.ones(field1.shape)*z[:,None,None])
        
        for traj in range(1,nt-1):
            
            divtraj[:,:,:,traj,traj] = div[:,:,:]
            # First calculate forward trajectories
            
            varout = ndimage.map_coordinates(u,[ztemp_f.ravel()*rdz,ytemp_f.ravel()*rdy,xtemp_f.ravel()*rdx],order=1,cval=np.nan)
            k1x = dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(v,[ztemp_f.ravel()*rdz,ytemp_f.ravel()*rdy,xtemp_f.ravel()*rdx],order=1,cval=np.nan)
            k1y = dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(w,[ztemp_f.ravel()*rdz,ytemp_f.ravel()*rdy,xtemp_f.ravel()*rdx],order=1,cval=np.nan)
            k1z = dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(u,[(ztemp_f+(0.5*k1z)).ravel()*rdz,(ytemp_f+(0.5*k1y)).ravel()*rdy,(xtemp_f+(0.5*k1x)).ravel()*rdx],order=1,cval=np.nan)
            k2x = dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(v,[(ztemp_f+(0.5*k1z)).ravel()*rdz,(ytemp_f+(0.5*k1y)).ravel()*rdy,(xtemp_f+(0.5*k1x)).ravel()*rdx],order=1,cval=np.nan)
            k2y = dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(w,[(ztemp_f+(0.5*k1z)).ravel()*rdz,(ytemp_f+(0.5*k1y)).ravel()*rdy,(xtemp_f+(0.5*k1x)).ravel()*rdx],order=1,cval=np.nan)
            k2z = dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(u,[(ztemp_f+(0.5*k2z)).ravel()*rdz,(ytemp_f+(0.5*k2y)).ravel()*rdy,(xtemp_f+(0.5*k2x)).ravel()*rdx],order=1,cval=np.nan)
            k3x = dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(v,[(ztemp_f+(0.5*k2z)).ravel()*rdz,(ytemp_f+(0.5*k2y)).ravel()*rdy,(xtemp_f+(0.5*k2x)).ravel()*rdx],order=1,cval=np.nan)
            k3y = dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(w,[(ztemp_f+(0.5*k2z)).ravel()*rdz,(ytemp_f+(0.5*k2y)).ravel()*rdy,(xtemp_f+(0.5*k2x)).ravel()*rdx],order=1,cval=np.nan)
            k3z = dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(u,[(ztemp_f+k3z).ravel()*rdz,(ytemp_f+k3y).ravel()*rdy,(xtemp_f+k3x).ravel()*rdx],order=1,cval=np.nan)
            k4x = dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(v,[(ztemp_f+k3z).ravel()*rdz,(ytemp_f+k3y).ravel()*rdy,(xtemp_f+k3x).ravel()*rdx],order=1,cval=np.nan)
            k4y = dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(w,[(ztemp_f+k3z).ravel()*rdz,(ytemp_f+k3y).ravel()*rdy,(xtemp_f+k3x).ravel()*rdx],order=1,cval=np.nan)
            k4z = dt * np.copy(varout.reshape(field1.shape))
            
            xtemp_f = xtemp_f + (k1x + 2.*k2x + 2.*k3x + k4x)/6.
            ytemp_f = ytemp_f + (k1y + 2.*k2y + 2.*k3y + k4y)/6.
            ztemp_f = ztemp_f + (k1z + 2.*k2z + 2.*k3z + k4z)/6.
            
            ind = np.arange(1,nt-traj)
            varout = ndimage.map_coordinates(div,[ztemp_f.ravel()*rdz,ytemp_f.ravel()*rdy,xtemp_f.ravel()*rdx],order=1,cval=np.nan)
            divtraj[:,:,:,ind,ind+traj] = np.array([(np.copy(varout.reshape(field1.shape))).T]*(nt-1-traj)).T
            
            varout = ndimage.map_coordinates(field2,[ztemp_f.ravel()*rdz,ytemp_f.ravel()*rdy,xtemp_f.ravel()*rdx],order=1,cval=np.nan)
            forwardref[:,:,:,nt-1-traj] = np.copy(varout.reshape(field1.shape))
            
            # Now do backward trajectories
            
            varout = ndimage.map_coordinates(u,[ztemp_b.ravel()*rdz,ytemp_b.ravel()*rdy,xtemp_b.ravel()*rdx],order=1,cval=np.nan)
            k1x = -dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(v,[ztemp_b.ravel()*rdz,ytemp_b.ravel()*rdy,xtemp_b.ravel()*rdx],order=1,cval=np.nan)
            k1y = -dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(w,[ztemp_b.ravel()*rdz,ytemp_b.ravel()*rdy,xtemp_b.ravel()*rdx],order=1,cval=np.nan)
            k1z = -dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(u,[(ztemp_b+(0.5*k1z)).ravel()*rdz,(ytemp_b+(0.5*k1y)).ravel()*rdy,(xtemp_b+(0.5*k1x)).ravel()*rdx],order=1,cval=np.nan)
            k2x = -dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(v,[(ztemp_b+(0.5*k1z)).ravel()*rdz,(ytemp_b+(0.5*k1y)).ravel()*rdy,(xtemp_b+(0.5*k1x)).ravel()*rdx],order=1,cval=np.nan)
            k2y = -dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(w,[(ztemp_b+(0.5*k1z)).ravel()*rdz,(ytemp_b+(0.5*k1y)).ravel()*rdy,(xtemp_b+(0.5*k1x)).ravel()*rdx],order=1,cval=np.nan)
            k2z = -dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(u,[(ztemp_b+(0.5*k2z)).ravel()*rdz,(ytemp_b+(0.5*k2y)).ravel()*rdy,(xtemp_b+(0.5*k2x)).ravel()*rdx],order=1,cval=np.nan)
            k3x = -dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(v,[(ztemp_b+(0.5*k2z)).ravel()*rdz,(ytemp_b+(0.5*k2y)).ravel()*rdy,(xtemp_b+(0.5*k2x)).ravel()*rdx],order=1,cval=np.nan)
            k3y = -dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(w,[(ztemp_b+(0.5*k2z)).ravel()*rdz,(ytemp_b+(0.5*k2y)).ravel()*rdy,(xtemp_b+(0.5*k2x)).ravel()*rdx],order=1,cval=np.nan)
            k3z = -dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(u,[(ztemp_b+k3z).ravel()*rdz,(ytemp_b+k3y).ravel()*rdy,(xtemp_b+k3x).ravel()*rdx],order=1,cval=np.nan)
            k4x = -dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(v,[(ztemp_b+k3z).ravel()*rdz,(ytemp_b+k3y).ravel()*rdy,(xtemp_b+k3x).ravel()*rdx],order=1,cval=np.nan)
            k4y = -dt * np.copy(varout.reshape(field1.shape))
            
            varout = ndimage.map_coordinates(w,[(ztemp_b+k3z).ravel()*rdz,(ytemp_b+k3y).ravel()*rdy,(xtemp_b+k3x).ravel()*rdx],order=1,cval=np.nan)
            k4z = -dt * np.copy(varout.reshape(field1.shape))
            
            xtemp_b = xtemp_b + (k1x + 2.*k2x + 2.*k3x + k4x)/6.
            ytemp_b = ytemp_b + (k1y + 2.*k2y + 2.*k3y + k4y)/6.
            ztemp_b = ztemp_b + (k1z + 2.*k2z + 2.*k3z + k4z)/6.
            
            ind = np.arange(traj,nt-1)
            varout = ndimage.map_coordinates(div,[ztemp_b.ravel()*rdz,ytemp_b.ravel()*rdy,xtemp_b.ravel()*rdx],order=1,cval=np.nan)
            divtraj[:,:,:,ind,ind-traj] = np.array([(np.copy(varout.reshape(field1.shape))).T]*(nt-1-traj)).T
            
            varout = ndimage.map_coordinates(field1,[ztemp_b.ravel()*rdz,ytemp_b.ravel()*rdy,xtemp_b.ravel()*rdx],order=1,cval=np.nan)
            backwardref[:,:,:,traj] = np.copy(varout.reshape(field1.shape))
            
        # Now you have to calculate the integral of the divergence of u and v
        # along the trajectories
        
        gral1 = np.zeros(divtraj.shape)
        
        gral1[:,:,:,:,1:] = scipy.integrate.cumtrapz(divtraj,dx=dt,axis=4)

        # and another intergral which is the integral of the exponential of the 
        # previous integral

        gral2 = np.zeros(divtraj.shape)

        gral2[:,:,:,:,1:] = scipy.integrate.cumtrapz(np.exp(-gral1),dx=dt,axis=4)
            
        # Here we construct the reflectivity field using the trajectories
        # and the integrals we just calculated
        
        ratio = gral2[:,:,:,np.arange(nt),np.arange(nt)]/gral2[:,:,:,np.arange(nt),-1]
        
        ref = backwardref + (ratio * (forwardref - backwardref))
        
        ref[:,:,:,0] = backwardref[:,:,:,0]
        ref[:,:,:,-1] = backwardref[:,:,:,-1]
            
        dRdt = 0.5 * rdt * (ref[1:ref.shape[0]-1,1:ref.shape[1]-1,1:ref.shape[2]-1,2:] - 
                            ref[1:ref.shape[0]-1,1:ref.shape[1]-1,1:ref.shape[2]-1,:ref.shape[3]-2])
        
        dRdz = 0.5 * rdz * (ref[2:,1:ref.shape[1]-1,1:ref.shape[2]-1,1:ref.shape[3]-1] -
                            ref[:ref.shape[0]-2,1:ref.shape[1]-1,1:ref.shape[2]-1,1:ref.shape[3]-1])
        
        dRdy = 0.5 * rdy * (ref[1:ref.shape[0]-1,2:,1:ref.shape[2]-1,1:ref.shape[3]-1] -
                            ref[1:ref.shape[0]-1,:ref.shape[1]-2,1:ref.shape[2]-1,1:ref.shape[3]-1])
        
        dRdx = 0.5 * rdx * (ref[1:ref.shape[0]-1,1:ref.shape[1]-1,2:,1:ref.shape[3]-1] -
                            ref[1:ref.shape[0]-1,1:ref.shape[1]-1,:ref.shape[2]-2,1:ref.shape[3]-1])
        
        foo = np.where((np.isnan(dRdt)) | (np.isnan(dRdy)) | (np.isnan(dRdx)) | (np.isnan(dRdz)))
        
        dRdt[foo] = np.nan
        dRdx[foo] = np.nan
        dRdy[foo] = np.nan
        dRdz[foo] = np.nan
        
        a = np.zeros(u.shape)
        b = np.zeros(u.shape)
        c = np.zeros(u.shape)
        d = np.zeros(u.shape)
        dw = np.zeros(u.shape)
        e = np.zeros(u.shape)
        f = np.zeros(u.shape)
        g = np.zeros(u.shape)
        gw = np.zeros(u.shape)
        h = np.zeros(u.shape)
        l = np.zeros(u.shape)
        
         
        a[1:ref.shape[0]-1,1:ref.shape[1]-1,1:ref.shape[2]-1] = np.nansum(dx * dx * dRdt * dRdx / (2.*(2.*beta+gamma)), axis = 3)/(nt-2)
        b[1:ref.shape[0]-1,1:ref.shape[1]-1,1:ref.shape[2]-1] = np.nansum(dx * dx * dRdx * dRdx / (2.*(2.*beta+gamma)), axis = 3)/(nt-2)
        c[1:ref.shape[0]-1,1:ref.shape[1]-1,1:ref.shape[2]-1] = np.nansum(dx * dx * dRdx * dRdy / (2.*(2.*beta+gamma)), axis = 3)/(nt-2)
        d[1:ref.shape[0]-1,1:ref.shape[1]-1,1:ref.shape[2]-1] = np.nansum(dx * dx * dRdx * dRdz / (2.*(2.*beta+gamma)), axis = 3)/(nt-2)
        dw[1:ref.shape[0]-1,1:ref.shape[1]-1,1:ref.shape[2]-1] = np.nansum(dx * dx * dRdx * dRdz / (2.*(2.*eta+nu)), axis = 3)/(nt-2)
        e[1:ref.shape[0]-1,1:ref.shape[1]-1,1:ref.shape[2]-1] = np.nansum(dx * dx * dRdt * dRdy / (2.*(2.*beta+gamma)), axis = 3)/(nt-2)    
        f[1:ref.shape[0]-1,1:ref.shape[1]-1,1:ref.shape[2]-1] = np.nansum(dx * dx * dRdy * dRdy / (2.*(2.*beta+gamma)), axis = 3)/(nt-2)
        g[1:ref.shape[0]-1,1:ref.shape[1]-1,1:ref.shape[2]-1] = np.nansum(dx * dx * dRdy * dRdz / (2.*(2.*beta+gamma)), axis = 3)/(nt-2)
        gw[1:ref.shape[0]-1,1:ref.shape[1]-1,1:ref.shape[2]-1] = np.nansum(dx * dx * dRdy * dRdz / (2.*(2.*eta+nu)), axis = 3)/(nt-2)
        h[1:ref.shape[0]-1,1:ref.shape[1]-1,1:ref.shape[2]-1] = np.nansum(dx * dx * dRdt * dRdz / (2.*(2.*eta+nu)), axis = 3)/(nt-2)
        l[1:ref.shape[0]-1,1:ref.shape[1]-1,1:ref.shape[2]-1] = np.nansum(dx * dx * dRdz * dRdz / (2.*(2.*eta+nu)), axis = 3)/(nt-2)
    
        # Now we call the function that performs the successive overrelaxation
        # to solve the PDEs for U and V. I am unaware of any Python packages
        # that have PDE solvers that would work here so this function is a 
        # f2py'ed fortran function
        
        uold = np.copy(u)
        vold = np.copy(v)
        wold = np.copy(w)
        
        u_new, v_new, w_new= relax3d(a,b,c,d,dw,e,f,g,gw,h,l,
                               u,v,w,beta,gamma,eta,nu,
                               itermax,tol,relax)
        
        # Check to see if u, v have converged in the big iteration
        
        dif = np.sqrt((u_new-uold)**2 + (v_new-vold)**2 + (w_new-wold)**2)
        if np.max(dif) < tol2:
            if verbose:
                print("*** Convergence in big iteration loop ***")
            break
        else:
            if verbose:
                print(" ")
                print("itermain = " + str(itermain))
                print("U relaxation solution in center of domain: " + str(u_new[int(u.shape[0]/2),int(u.shape[1]/2),int(u.shape[2]/2)]))
                print("V relaxation solution in center of domain: " + str(v_new[int(u.shape[0]/2),int(u.shape[1]/2),int(u.shape[2]/2)]))
                print("W relaxation solution in center of domain: " + str(w_new[int(u.shape[0]/2),int(u.shape[1]/2),int(u.shape[2]/2)]))
                print("Maximum difference: " + str(np.max(dif)))
        
            u = (1.0 - under) * uold + under * u_new
            v = (1.0 - under) * vold + under * v_new
            w = (1.0 - under) * wold + under * w_new
            
            div[1:div.shape[0]-1,1:div.shape[1]-1,1:div.shape[2]-1] = (
                0.5*rdz*(w[2:,1:w.shape[1]-1,1:w.shape[2]-1] - w[:w.shape[0]-2,1:w.shape[1]-1,1:w.shape[2]-1]) +
                0.5*rdy*(v[1:v.shape[0]-1,2:,1:v.shape[2]-1] - v[1:v.shape[0]-1,:v.shape[1]-2,1:v.shape[1]-1]) +
                0.5*rdx*(u[1:u.shape[0]-1,1:u.shape[1]-1,2:] - u[1:u.shape[0]-1,1:u.shape[1]-1,:u.shape[2]-2])
                )
    
            div[:,0,:] = div[:,1,:]
            div[:,-1,:] = div[:,-2,:]
    
            div[:,:,0] = div[:,:,1]
            div[:,:,-1] = div[:,:,-2]
    
            div[0,:,:] = div[1,:,:]
            div[-1,:,:] = div[-2,:,:]
        
    return u, v, w, ref