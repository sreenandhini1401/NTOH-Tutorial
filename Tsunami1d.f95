
program Tsunami
implicit none

real,parameter:: g=9.81
integer::x,t,n,i,j
real::del_x,del_t,a,c,dt
real,allocatable,dimension(:,:):: eta,vel,m,h
!real,dimension(801)::x1



x=100
del_x=0.125
!x1= [0, (i*1/8 ,i=1,799),x]
t=3000
n=(x/del_x)+1
print*,n

allocate(eta(t,n),m(t,n),h(t,n))

eta=0
m=0
h=0

do j=1,t
do i=1,n
	h(j,i)= int(50-45*tanh((del_x*(i-1)-70)/8))
end do
end do

!dt = del_x/sqrt(2*9.81*95);
del_t=0.001;
a= del_t/(2*del_x);
c= del_t/del_x;

!!!!Initial condition
do i=2,n-1
    eta(1,i)=0.5*exp(-((del_x*(i-1)-20)**2)/8)
end do
	eta(1,1)=eta(1,2)
	eta(1,n)=eta(1,n-1)
do i=2,n-1
    m(1,i)=100*eta(1,i)
end do
	m(1,1)=m(1,2)
	m(1,n)=m(1,n-1)
!!!!!FTCS
        do j= 2,n-1
            m(2,j)= m(1,j) - a*g*(eta(1,j)+h(2,j))*(eta(1,j+1)-eta(1,j-1));
        end do
        m(2,1)=0
		m(2,n)=0
		do j= 2,n-1
            eta(2,j)= eta(1,j)-a*(m(1,j+1)-m(1,j-1));
        end do
		eta(2,1)=eta(2,2)
		eta(2,n)=eta(2,n-1)

!!!!!!!!!leap frog
		do i=3,t
			do j=2,n-1

				m(i,j)= m(i-2,j)- g*c*((eta(i-1,j+1)+eta(i-1,j))+h(i-1,j))*(eta(i-1,j+1)-eta(i-1,j-1))
			end do
			m(i,1)=m(i,2)
			m(1,n)=m(1,n-1)
			do j=2,n-1
				eta(i,j)= eta(i-2,j)- c*(m(i-1,j+1)-m(i-1,j-1));

			end do
			eta(i,1)=eta(i,2)
			eta(i,n)=eta(i,n-1)

		end do

		open(20,file='eta_1d.txt',status='unknown')
		do i=1,n
			write(20,*) eta(:,i)
		end do
		close(20)

		open(20,file='m_1d.txt',status='unknown')
		do i=1,n
			write(20,*) m(:,i)
		end do
		close(20)
end

