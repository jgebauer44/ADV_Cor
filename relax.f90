      subroutine relax2d(a,b,c,d,e,u,v,ny,nx,itermax,tol, &
                         relax)
         implicit none

         integer, intent(in)  :: nx, ny, itermax
         real,   intent(in)  :: a(ny,nx), b(ny,nx), c(ny,nx)
         real,   intent(in)  :: d(ny,nx), e(ny,nx)
         real,   intent(in)  :: tol
         real,   intent(inout) :: u(ny,nx), v(ny,nx)
!f2py intent(in,out) :: u, v
         real,   intent(in) :: relax
         real uold(ny,nx), vold(ny,nx)
         real utemp, vtemp, maxdif, ubar, vbar
         real dif
         integer icount, iter, i, j

         do j = 1, nx
            do i = 1, ny
               uold(i,j) = u(i,j)
               vold(i,j) = v(i,j)
            end do
         end do

         icount = 0
         do iter = 1, itermax
            icount = icount + 1
            
            if (icount/10 .eq. 1) then
               maxdif = 0.
               do j = 1, nx
                  do i = 1, ny
                     dif = sqrt((u(i,j) - uold(i,j))**2 + &
                                (v(i,j) - vold(i,j))**2)
                     maxdif = max(maxdif,dif)
                     uold(i,j) = u(i,j)
                     vold(i,j) = v(i,j)
                   end do
                end do
 
                icount = 0

                if (maxdif .lt. tol) exit

            end if
        
          do j = 2, nx-1
           do i = 2, ny-1

            utemp = u(i,j)
            vtemp = v(i,j)

            ubar = 0.25 * (u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1))
            vbar = 0.25 * (v(i+1,j) + v(i-1,j) + v(i,j+1) + v(i,j-1))

            u(i,j) = ubar - 0.25*(a(i,j) + c(i,j)*v(i,j))
            u(i,j) = u(i,j)/(1. + 0.25*b(i,j))

            u(i,j) = (1.0-relax)*utemp + relax*u(i,j)
          
            v(i,j) = vbar - 0.25*(d(i,j) + c(i,j)*u(i,j) )
            v(i,j) = v(i,j)/(1. + 0.25*e(i,j))
 
            v(i,j) = (1.0 - relax)*vtemp + relax*v(i,j)
           end do
          end do

          do i = 2, ny-1
            u(i,1) = u(i,2)
            u(i,nx) = u(i,nx-1)

            v(i,1) = v(i,2)
            v(i,nx) = v(i,nx-1)
          end do

          do j = 1, nx
            u(1,j) = u(2,j)
            u(ny,j) = u(ny-1,j)

            v(1,j) = v(2,j)
            v(ny,j) = v(ny-1,j)
          end do
         end do
      end subroutine
 
           

      subroutine relax3d(a,b,c,d,dw,e,f,g,gw,h,l,u,v,w, &
                        nz,ny,nx,beta,gamma,eta,nu, &
                        itermax,tol,relax)
         implicit none

         integer, intent(in)  :: nx, ny,nz, itermax
         real,   intent(in)  :: a(nz,ny,nx), b(nz,ny,nx) 
         real,   intent(in)  :: c(nz,ny,nx), d(nz,ny,nx)
         real,   intent(in)  :: dw(nz,ny,nx), e(nz,ny,nx)
         real,   intent(in)  :: f(nz,ny,nx), g(nz,ny,nx)
         real,   intent(in)  :: gw(nz,ny,nx), h(nz,ny,nx)
         real,   intent(in)  :: l(nz,ny,nx)
         real,   intent(in)  :: tol
         real,   intent(in)  :: relax
         real,   intent(in)  :: beta,gamma,eta,nu
         real,   intent(inout) :: u(nz,ny,nx), v(nz,ny,nx)
         real,   intent(inout)  ::  w(nz,ny,nx)
!f2py intent(in,out) :: u, v, w
         real uold(nz,ny,nx), vold(nz,ny,nx), wold(nz,ny,nx)
         real utemp, vtemp, wtemp, maxdif, ubar1, vbar1, wbar1
         real ubar2, vbar2, wbar2, dif
         integer icount, iter, i, j, k

         do k = 1, nx
          do j = 1, ny
            do i = 1, nz
               uold(i,j,k) = u(i,j,k)
               vold(i,j,k) = v(i,j,k)
               wold(i,j,k) = w(i,j,k)
            end do
           end do
         end do

         icount = 0
         do iter = 1, itermax
            icount = icount + 1
            
            if (icount/10 .eq. 1) then
               maxdif = 0.
               do k = 1, nx
                do j = 1, ny
                  do i = 1, nz
                     dif = sqrt((u(i,j,k) - uold(i,j,k))**2 + &
                                (v(i,j,k) - vold(i,j,k))**2 + &
                                (w(i,j,k) - wold(i,j,k))**2)
                     maxdif = max(maxdif,dif)
                     uold(i,j,k) = u(i,j,k)
                     vold(i,j,k) = v(i,j,k)
                     wold(i,j,k) = w(i,j,k)
                   end do
                 end do
                end do
 
                icount = 0

                if (maxdif .lt. tol) exit

            end if
        
          do k = 2, nx-1
           do j = 2, ny-1
            do i = 2, nz-1

             utemp = u(i,j,k)
             vtemp = v(i,j,k)
             wtemp = w(i,j,k)

             ubar1 = beta*(u(i,j,k+1) + u(i,j,k-1) + &
                        u(i,j+1,k) + u(i,j-1,k))/(2.*(2.*beta+gamma))

             ubar2 = gamma*(u(i+1,j,k) + u(i-1,j,k)) / &
                    (2*(2.*beta+gamma))

             vbar1 = beta*(v(i,j,k+1) + v(i,j,k-1) + &
                        v(i,j+1,k) + v(i,j-1,k))/(2.*(2.*beta+gamma))

             vbar2 = gamma*(v(i+1,j,k) + v(i-1,j,k)) / &
                    (2*(2.*beta+gamma))

             wbar1 = eta*(w(i,j,k+1) + w(i,j,k-1) + &
                        w(i,j+1,k) + w(i,j-1,k))/(2.*(2.*eta+nu))

             wbar2 = nu*(w(i+1,j,k) + w(i-1,j,k)) / &
                    (2*(2.*eta+nu))

             u(i,j,k) = ubar1 + ubar2 - (a(i,j,k) +  &
                                      c(i,j,k) * v(i,j,k) + &
                                      d(i,j,k) * w(i,j,k))

             u(i,j,k) = u(i,j,k)/(1. + b(i,j,k))

             u(i,j,k) = (1.0-relax)*utemp + relax*u(i,j,k)
            
             v(i,j,k) = vbar1 + vbar2 - (e(i,j,k) +  &
                                      c(i,j,k) * u(i,j,k) + &
                                      g(i,j,k) * w(i,j,k))

             v(i,j,k) = v(i,j,k)/(1. + f(i,j,k))

             v(i,j,k) = (1.0 - relax)*vtemp + relax*v(i,j,k)
            
             w(i,j,k) = wbar1 + wbar2 - (h(i,j,k) +  &
                                      dw(i,j,k) * u(i,j,k) + &
                                      gw(i,j,k) * v(i,j,k))

             w(i,j,k) = w(i,j,k)/(1. + l(i,j,k))

             w(i,j,k) = (1.0-relax)*wtemp + relax*w(i,j,k)

            end do
           end do
          end do

          do i = 2, nz-1
           do k = 2, nx-1 
             u(i,1,k) = u(i,2,k)
             u(i,ny,k) = u(i,ny-1,k)

             v(i,1,k) = v(i,2,k)
             v(i,ny,k) = v(i,ny-1,k)

             w(i,1,k) = w(i,2,k)
             w(i,ny,k) = w(i,ny-1,k)
           end do
          end do

          do j = 1, ny
           do i = 2,nz-1
            u(i,j,1) = u(i,j,2)
            u(i,j,nx) = u(i,j,nx-1)
            
            v(i,j,1) = v(i,j,2)
            v(i,j,nx) = v(i,j,nx-1)

            w(i,j,1) = w(i,j,2)
            w(i,j,nx) = w(i,j,nx-1)
           end do
          end do

          do k = 1, nx
           do j = 1, ny
            u(1,j,k) = u(2,j,k)
            u(nz,j,k) = u(nz-1,j,k)

            v(1,j,k) = v(2,j,k)
            v(nz,j,k) = v(nz-1,j,k)

            w(1,j,k) = w(2,j,k)
            w(nz,j,k) = w(nz-1,j,k)
           end do
          end do
              
         end do
      end subroutine

