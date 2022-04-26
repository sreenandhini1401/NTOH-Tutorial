subroutine sinft2d(data,n1,n2,nmax,isign)
    implicit none
      integer:: n1,n2,nmax,isign
      real:: data(n1,n2)
      integer::i,j
      real:: wk(nmax),temp

      do  i=1,n1
        do  j=1,n2
            wk(j)=data(i,j)
        end do
        call sinft(wk,n2)
        do  j=1,n2
            data(i,j)=wk(j)
        end do
      end do

      do  j=1,n2
        do  i=1,n1
            wk(i)=data(i,j)
        end do
        call sinft(wk,n1)
        do  i=1,n1
            data(i,j)=wk(i)
        end do
      end do


      if (isign.lt.0) then
      temp=4.0d0/n1/n2
        do  i=1,n1
            do  j=1,n2
                data(i,j)=data(i,j)*temp
            end do
        end do
      end if
      return
      end


subroutine sinft(y,n)
    implicit none
      integer:: n
      real:: y(n)
      integer:: j
      real:: sum,y1,y2
      real:: theta,wi,wpi,wpr, wr,wtemp
!     Double precision in the trigonometric recurrences.
      theta=3.141592653589793d0/dble(n)   ! Initialize the recurrence.
      wr=1.0d0
      wi=0.0d0
      wpr=-2.0d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      y(1)=0.0d0
    do  j=1,n/2
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
          y1=wi*(y(j+1)+y(n-j+1))
          y2=0.5d0*(y(j+1)-y(n-j+1))
          y(j+1)=y1+y2
          y(n-j+1)=y1-y2
    end do
      call realft(y,n,1)
      sum=0.0d0
      y(1)=y(1)/2.0d0
      y(2)=0.0d0
    do  j=1,n-1,2
      sum=sum+y(j)
      y(j)=y(j+1)
      y(j+1)=sum
    end do
      return
      end



      subroutine realft(data,n,isign)
        implicit none
      integer:: isign,n
      real:: data(n)
! USES four1 Calculates the Fourier transform of a set
! of n real-valued data points. Replaces this data
!(which is stored in array data(1:n)) by the positive
! frequency half of its complex Fourier transform. The
! real-valued  rst and last components of the complex
! transform are returned as elements data(1) and data(2),
! respectively. n must be a power of 2. This routine also
! calculates the inverse transform of a complex data array
! if it is the transform of real data. (Result in this case
! must be multiplied by 2/n.)
      INTEGER:: i,i1,i2,i3,i4,n2p3
      real:: c1,c2,h1i,h1r,h2i,h2r,wis,wrs
      real:: theta,wi,wpi,wpr, wr,wtemp
!      Double precision for the trigonometric recurrences.
      theta=3.141592653589793d0/dble(n/2) !Initialize the recurrence.
      c1=0.5d0
      if (isign.eq.1) then
      c2=-0.5d0
      call four1(data,n/2,1)   !The forward transform is here.
      else
      c2=0.5d0   !Otherwise set up for an inverse transform.
      theta=-theta
      end if
      wpr=-2.0d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      wr=1.0d0+wpr
      wi=wpi
      n2p3=n+3
      do  i=2,n/4  ! Case i=1 done separately below.
      i1=2*i-1
      i2=i1+1
      i3=n2p3-i2
      i4=i3+1
!      wrs=sngl(wr)
!      wis=sngl(wi)
      wrs=(wr)
      wis=(wi)
      h1r=c1*(data(i1)+data(i3))
! The two separate transforms are separated out of data.
      h1i=c1*(data(i2)-data(i4))
      h2r=-c2*(data(i2)+data(i4))
      h2i=c2*(data(i1)-data(i3))
      data(i1)=h1r+wrs*h2r-wis*h2i
! Here they are recombined to form the true trans- form of the original real data.
      data(i2)=h1i+wrs*h2i+wis*h2r
      data(i3)=h1r-wrs*h2r+wis*h2i
      data(i4)=-h1i+wrs*h2i+wis*h2r
      wtemp=wr            !The recurrence.
      wr=wr*wpr-wi*wpi+wr
      wi=wi*wpr+wtemp*wpi+wi
    end do
      if (isign.eq.1) then
      h1r=data(1)
      data(1)=h1r+data(2)
      data(2)=h1r-data(2)
! Squeeze the  rst and last data together to get them all within the original array.
      else
      h1r=data(1)
      data(1)=c1*(h1r+data(2))
      data(2)=c1*(h1r-data(2))
      call four1(data,n/2,-1)  !This is the inverse transformfor the case isign=-1.
      end if
      return
      END


      subroutine  four1(data,nn,isign)
        implicit none
      integer:: isign,nn
      real:: data(2*nn)
! Replaces data(1:2*nn) by its discrete Fourier transform,
! if isign is input as 1; or replaces data(1:2*nn) by nn
! times its inverse discrete Fourier transform, if isign
! is input as  -1. data is a complex array of length nn or,
! equivalently, a real array of length 2*nn. nn MUST be an
! integer power of 2 (this is not checked for!).
      integer:: i,istep,j,m,mmax,n
      real:: tempi,tempr
      real:: theta,wi,wpi,wpr,wr,wtemp
!      Double precision for the trigonometric recurrences.
      n=2*nn
      j=1
      do  i=1,n,2  !This is the bit-reversal section of the routine.
      if(j.gt.i)then
      tempr=data(j)  !Exchange the two complex numbers.
      tempi=data(j+1)
      data(j)=data(i)
      data(j+1)=data(i+1)
      data(i)=tempr
      data(i+1)=tempi
      end if
      m=n/2
   1  if ((m.ge.2).and.(j.gt.m))  then
      j=j-m
      m=m/2
      goto 1
      end if
      j=j+m
    end do
      mmax=2  !Here begins the Danielson-Lanczos section of the routine.
   2  if (n.gt.mmax) then  !Outer loop executed log 2 nn times.
      istep=2*mmax
      theta=6.28318530717959d0/(isign*mmax) !Initialize for the trigonometric recurrence.
      wpr=-2.d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      wr=1.d0
      wi=0.d0
      do  m=1,mmax,2  !Here are the two nested inner loops.
      do i=m,n,istep
      j=i+mmax  !This is the Danielson-Lanczos formula:
!      tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
!     tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
      tempr=(wr)*data(j)-(wi)*data(j+1)
      tempi=(wr)*data(j+1)+(wi)*data(j)
      data(j)=data(i)-tempr
      data(j+1)=data(i+1)-tempi
      data(i)=data(i)+tempr
      data(i+1)=data(i+1)+tempi
 end do
      wtemp=wr !Trigonometric recurrence.
      wr=wr*wpr-wi*wpi+wr
      wi=wi*wpr+wtemp*wpi+wi
 end do
      mmax=istep
      goto 2 !Not yet done.
      end if !All done.
      return
      END

    subroutine output(filename,un,n,variable)
    implicit none
    integer:: un,i
    character (len=50):: filename
    real::n
    real,dimension(int(n),int(n)):: variable

    open (un,file=filename,status='unknown')
    do i= 1,n
        write(un,*),variable(i,:)
    end do

end subroutine
