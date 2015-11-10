!-----------------------------------------------------------------------------!
! MAE 5093 DOCUMENTATION | Engineering Numerical Analysis
!-----------------------------------------------------------------------------!
! >>> Fast Poisson Solver
!     Fast Sine Solver for solving 2D Poisson equation
!     d2u/dx2 + d2u/dy2 = f(x,y)
!     Drichlet b.c.
!
!-----------------------------------------------------------------------------!
! References: 
! * Fundamentals of Engineering Numerical Analysis by P. Moin (2012)
! * Numerical Recipes: The Art of Scientific Computing, 2nd Edition (1992) 
!-----------------------------------------------------------------------------!
! Written by Omer San
!            CFDLab, Oklahoma State University, cfdlab.osu@gmail.com
!            www.cfdlab.org
! 
! Last updated: Nov. 11, 2015
!-----------------------------------------------------------------------------!

program poisson2d
implicit none
integer::i,j,nx,ny
real*8,dimension(:,:),allocatable ::u,f,ue,e
real*8,dimension(:),allocatable ::x,y
real*8 ::dx,dy,x0,xL,y0,yL

!Domain
x0 =-1.0d0 !left
xL = 1.0d0 !right

y0 =-1.0d0 !bottom
yL = 1.0d0 !up

!number of points 
nx = 64 !number of grid points in x (i.e., should be power of 2)
ny = nx  !number of grid points in y

!grid spacing (spatial)
dx = (xL-x0)/dfloat(nx)
dy = (yL-y0)/dfloat(ny)

!spatial coordinates 
allocate(x(0:nx))
do i=0,nx
x(i) = x0 + dfloat(i)*dx
end do

allocate(y(0:ny))
do j=0,ny
y(j) = y0 + dfloat(j)*dy
end do


allocate(u(0:nx,0:ny))
allocate(f(0:nx,0:ny))
allocate(e(0:nx,0:ny))
allocate(ue(0:nx,0:ny))

!---------------------------------------------!
!Exact solution (test case from Moin's textbook):
!---------------------------------------------!
do j=0,ny
do i=0,nx
f(i,j) =-2.0d0*(2.0d0-x(i)*x(i)-y(j)*y(j))
ue(i,j)= (x(i)*x(i)-1.0d0)*(y(j)*y(j)-1.0d0)
end do
end do


!Drichlet boundary conditions
do i=0,nx
u(i,0)  = 0.0d0	
u(i,ny) = 0.0d0				  					  	
end do

do j=0,ny
u(0,j)  = 0.0d0		
u(nx,j) = 0.0d0					  	
end do

!----------------------!
!Fast Poisson Solver:
!----------------------!
call FPS(nx,ny,dx,dy,f,u)

!----------------------!
!Error analysis:
!----------------------!
do i=0,nx
do j=0,ny
e(i,j) = dabs(u(i,j)-ue(i,j))
end do 
end do

!maximum norm
write(*,*)"Max-norm =",maxval(e)


!Plot field
open(10,file='field.plt')
write(10,*) 'variables ="x","y","f","u","ue"'
write(10,*)'zone f=point i=',nx+1,',j=',ny+1
do j=0,ny
do i=0,nx
write(10,*) x(i),y(j),f(i,j),u(i,j),ue(i,j)
end do
end do
close(10)

end


!---------------------------------------------------------------------------!
!Routines for fast Poisson solver
!---------------------------------------------------------------------------!

!---------------------------------------------------------------------------!
!fast poisson solver (FPS)
!
!fast direct poisson solver for homogeneous drichlet boundary conditions
!using discreate fast sin transformation along x and y axis 
!underlying second order finite difference formula
!---------------------------------------------------------------------------!
subroutine FPS(nx,ny,dx,dy,f,u) 
implicit none
integer::i,j,nx,ny,isign
real*8,dimension(0:nx,0:ny)::u,f
real*8,dimension(:,:),allocatable:: ft
real*8::dx,dy,pi,alpha

pi=4.0d0*datan(1.0d0)

allocate(ft(0:nx,0:ny))

!rename for rhs (not to replace the data)
do i=0,nx
do j=0,ny
ft(i,j) = f(i,j)
end do
end do

!fast inverse fourier sine transform of source term:
isign=-1
call sinft2d(nx,ny,isign,ft)

!Compute fourier coefficient of u:
!Homegenous Drichlet bc applies 
do i=1,nx-1
do j=1,ny-1
alpha=2.0d0/(dx*dx)*(dcos(pi*dfloat(i)/dfloat(nx))-1.0d0) &
     +2.0d0/(dy*dy)*(dcos(pi*dfloat(j)/dfloat(ny))-1.0d0)

u(i,j)=ft(i,j)/alpha
end do
end do

!fast forward fourier sine transform:
isign=1
call sinft2d(nx,ny,isign,u)

deallocate(ft)

return
end



!---------------------------------------------------------------------------!
!Compute fast fourier sine transform for 2D data
!Homogeneous Drichlet Boundary Conditios (zero all boundaries)
!Input:: u(0:nx,0:ny) 
!        where indices 0,nx,ny represent boundary data and should be zero
!Output::override
!Automatically normalized
!isign=-1 is inverse transform and 2/N is already applied 
!        (from grid data to fourier coefficient)
!isign=+1 is forward transform
!        (from fourier coefficient to grid data)
!---------------------------------------------------------------------------!
subroutine sinft2d(nx,ny,isign,u)
implicit none
integer::nx,ny,isign
real*8 ::u(0:nx,0:ny)
integer::i,j
real*8, dimension(:),allocatable  :: v

if (isign.eq.-1) then !inverse transform
! compute inverse sine transform to find fourier coefficients of f in x-direction
allocate(v(nx))
do j=1,ny-1
	do i=1,nx
	v(i) = u(i-1,j)
	end do
	call sinft(v,nx)
	do i=2,nx
	u(i-1,j)=v(i)*2.0d0/dfloat(nx)
	end do
end do
deallocate(v)
allocate(v(ny))
! compute inverse sine transform to find fourier coefficients of f in y-direction
do i=1,nx-1
	do j=1,ny
	v(j) = u(i,j-1)
	end do
	call sinft(v,ny)
	do j=2,ny
	u(i,j-1)=v(j)*2.0d0/dfloat(ny)
	end do
end do
deallocate(v)
else  !forward transform
! compute forward sine transform to find fourier coefficients of f in x-direction
allocate(v(nx))
do j=1,ny-1
	do i=1,nx
	v(i) = u(i-1,j)
	end do
	call sinft(v,nx)
	do i=2,nx
	u(i-1,j)=v(i)
	end do
end do
deallocate(v)
allocate(v(ny))
! compute forward sine transform to find fourier coefficients of f in y-direction
do i=1,nx-1
	do j=1,ny
	v(j) = u(i,j-1)
	end do
	call sinft(v,ny)
	do j=2,ny
	u(i,j-1)=v(j)
	end do
end do
deallocate(v)
end if

return
end 

!---------------------------------------------------------------------------!
!calculates sine transform of a set of n real valued data points, y(1,2,..n)
!y(1) is zero, y(n) not need to be zero, but y(n+1)=0
!also calculates inverse transform, but output should be multiplied by 2/n
!n should be powers of 2
!use four1 and realft routines
!---------------------------------------------------------------------------!
subroutine sinft(y,n)
implicit none
integer::n,j,m
real*8 ::wr,wi,wpr,wpi,wtemp,theta,y1,y2,sum
real*8 ::y(n)
      theta=3.14159265358979d0/dble(n)
      wr=1.0d0
      wi=0.0d0
      wpr=-2.0d0*dsin(0.5d0*theta)**2
      wpi=dsin(theta)
      y(1)=0.0
      m=n/2
      do j=1,m
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
        y1=wi*(y(j+1)+y(n-j+1))
        y2=0.5*(y(j+1)-y(n-j+1))
        y(j+1)=y1+y2
        y(n-j+1)=y1-y2
      end do
      
      call realft(y,m,+1)
      
      sum=0.0
      y(1)=0.5*y(1)
      y(2)=0.0
      
      do j=1,n-1,2
        sum=sum+y(j)
        y(j)=y(j+1)
        y(j+1)=sum
      end do
return
end

!---------------------------------------------------------------------------!
!computes real fft
!---------------------------------------------------------------------------!
subroutine realft(data,n,isign)
implicit none
integer::n,isign,i,i1,i2,i3,i4,n2p3
real*8 ::wr,wi,wpr,wpi,wtemp,theta,c1,c2,h2r,h2i,h1r,h1i
real*8 ::data(*)
real   ::wrs,wis 
      theta=6.28318530717959d0/2.0d0/dble(n)
      c1=0.5
      if (isign.eq.1) then
        c2=-0.5
        call four1(data,n,+1)
      else
        c2=0.5
        theta=-theta
      endif
      wpr=-2.0d0*dsin(0.5d0*theta)**2
      wpi=dsin(theta)
      wr=1.0d0+wpr
      wi=wpi
      n2p3=2*n+3
      do i=2,n/2+1
        i1=2*i-1
        i2=i1+1
        i3=n2p3-i2
        i4=i3+1
        wrs=sngl(wr)
        wis=sngl(wi)
        h1r=c1*(data(i1)+data(i3))
        h1i=c1*(data(i2)-data(i4))
        h2r=-c2*(data(i2)+data(i4))
        h2i=c2*(data(i1)-data(i3))
        data(i1)=h1r+wrs*h2r-wis*h2i
        data(i2)=h1i+wrs*h2i+wis*h2r
        data(i3)=h1r-wrs*h2r+wis*h2i
        data(i4)=-h1i+wrs*h2i+wis*h2r
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
	  end do
      if (isign.eq.1) then
        h1r=data(1)
        data(1)=h1r+data(2)
        data(2)=h1r-data(2)
      else
        h1r=data(1)
        data(1)=c1*(h1r+data(2))
        data(2)=c1*(h1r-data(2))
        call four1(data,n,-1)
      endif
return
end

!---------------------------------------------------------------------------!
!FFT routine for 1-dimensional data 
!---------------------------------------------------------------------------!
subroutine four1(data,nn,isign)
implicit none
integer:: nn,isign,i,j,m,n,mmax,istep
real*8 :: wr,wi,wpr,wpi,wtemp,theta,tempr,tempi
real*8 :: data(*)
      n=2*nn
      j=1
      do i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        go to 1
        endif
        j=j+m
     end do
      mmax=2
2       if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*dsin(0.5d0*theta)**2
        wpi=dsin(theta)
        wr=1.d0
        wi=0.d0
        do m=1,mmax,2
          do i=m,n,istep
            j=i+mmax
            tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
            tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
		  end do
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
		end do
        mmax=istep
      go to 2
      endif
return
end



