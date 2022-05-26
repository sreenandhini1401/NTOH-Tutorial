program v_str_function

    implicit none
    integer::i,j,k,g,l,it,iii,n,caseno
    real:: length, h,dt,vel,n_iter,iter,error,t
    real:: max_res,beta,er,residual,re,tol,err,convergence,tol_1
    real:: conv_sf, c_v,c_sf,cfl,gamma1
    real::a,b,c,d
    real,allocatable:: si(:,:),zeta(:,:),P(:,:),u(:,:),v(:,:),si0(:,:),xg(:,:),yg(:,:)
    real,allocatable:: si_old(:,:),zeta_old(:,:),z_err(:,:),rhs(:,:),v1(:,:)

    print*,'CASE1: N=51,RE=100;CASE2: N=51,RE=400; CASE3: N=51,RE=1000 ;CASE4: N=101,RE=100'
    print*, 'Desired case: Select between 1 to 4 '
    read*,caseno
    print*,'Enter Length :'
    read*, length
    print*,'Enter Lid Cap Velocity :'
    read*, vel
    if ((caseno== 1) .or. (caseno ==2) .or. (caseno==3)) then
            n=51
    else
            n=101
    end if

    if ((caseno== 1) .or. (caseno ==4)) then
        re=100
    elseif (caseno==2) then
        re=400
    elseif (caseno==3) then
        re=1000
    end if

    h= length/(n-1)
    dt=0.001
	n_iter=10000
	t=0
    cfl= vel*dt/h
    gamma1= (dt)/(re*h**2)
	print*,cfl,gamma1

	error=1e-5
    tol=0.1
    convergence=0
    it=0
    beta=1.5
    allocate(si(n,n),zeta(n,n),p(n,n),u(n,n),v(n,n),xg(n,n),yg(n,n))
    allocate(si_old(n,n),zeta_old(n,n),z_err(n,n),rhs(n,n))
    si=0
    zeta=0
    P=0
    v=0
    u=0
    u(1,:)=vel
    rhs=0

    do iii=1,100000

		do iter=1,n_iter
            si_old= si
			do i=2,n-1
				 do j=2,n-1
					 si(i,j)= beta*0.25*(si(i+1,j)+si(i-1,j)+si(i,j+1)+&
					 si(i,j-1)+zeta(i,j)*h**2)+(1-beta)*si(i,j)
				end do
			end do

			call residue(si,si_old,n,h,c_sf)
			if (c_sf .le. error) then
                exit
            end if

		end do

    !!!Boundary conditions
    do i= 2,n-1
        do j= 2,n-1
            zeta(1,j)= -2.0*(si(2,j)/h**2)-(((2.0*vel))/h)
            zeta(n,j)= -2.0*(si(n-1,j))/h**2
            zeta(i,1)= -2.0*(si(i,2))/h**2
            zeta(i,n)= -2.0*(si(i,n-1))/h**2
        end do
    end do

    !!!!!vorticity equation
    do i= 2,n-1
        do j=2,n-1
            zeta(i,j)= zeta(i,j)- (0.25*dt/h**2)*(si(i+1,j)-si(i-1,j))*&
            (zeta(i,j+1)-zeta(i,j-1))+&
            (0.25*dt/h**2)*(si(i,j+1)-si(i,j-1))*(zeta(i+1,j)-zeta(i-1,j))&
            +(dt/(re*(h**2)))*(zeta(i+1,j)+zeta(i-1,j)+&
            zeta(i,j+1)+zeta(i,j-1)-4*zeta(i,j))
        end do
    end do

    u(1,:)= vel
    do i= 2,n-1
        do j= 2,n-1
            u(i,j)= -0.5*(si(i+1,j)-si(i-1,j))/(h)
            v(i,j)= (0.5*(si(i,j+1)-si(i,j-1)))/(h)
        end do
    end do
!!!pressure_bc
    p(n,1)=1.0
    do j=2,n
        a= -(1/(re*2*h))*(-3*zeta(n,j)+4*zeta(n-1,j)-zeta(n-2,j))
        p(n,j)= p(n,j-1)+a*h  !!forward difference
    end do
    do i=n-1,1,-1
        b= -(1/(re*2*h))*(-3*zeta(i,n)+4*zeta(i,n-1)-zeta(i,n-2))
        p(i,n)= p(i+1,n)+b*h
    end do
    do j=n-1,1,-1
        c= -(1/(re*2*h))*(-3*zeta(1,j)+4*zeta(2,j)-zeta(3,j))
        p(1,j)=p(1,j+1)+c*h
    end do
    do i=2,n-1
        d = -(1/(re*2*h))*(-3*zeta(i,1)+4*zeta(i,2)-zeta(i,3))
        p(i,1)= p(i-1,1)+ d*h
    end do
    !!pressure
    do i= 2,n-1
        do j=2,n-1
            rhs(i,j)= (((si(i,j+1)-2*si(i,j)+si(i,j-1))*(si(i+1,j)-2*si(i,j)+si(i-1,j))/(h**4))&
            -(((si(i+1,j+1)-si(i+1,j-1)-si(i-1,j+1)+si(i-1,j-1))/(4*h**2))**2))
            p(i,j)= 0.25*(p(i+1,j)+p(i-1,j)+p(i,j+1)+p(i,j-1))- 0.5*(rhs(i,j))
        end do
    end do
    it= it+1
	t=t+dt
end do

    if (caseno.eq.1) then

            call output('z1.csv',10,n,zeta)
            call output('s1.csv',30,n,si)
            call output('u1.csv',40,n,u)
            call output('v1.csv',50,n,v)
            call output('P1.csv',60,n,p)

    elseif (caseno.eq.2) then

            call output('z2.csv',10,n,zeta)
            call output('s2.csv',30,n,si)
            call output('u2.csv',40,n,u)
            call output('v2.csv',50,n,v)
            call output('P2.csv',60,n,p)

    elseif (caseno.eq.3) then

            call output('z3.csv',10,n,zeta)
            call output('s3.csv',30,n,si)
            call output('u3.csv',40,n,u)
            call output('v3.csv',50,n,v)
            call output('P3.csv',60,n,p)

     elseif (caseno.eq.4) then

            call output('z4.csv',10,n,zeta)
            call output('s4.csv',30,n,si)
            call output('u4.csv',40,n,u)
            call output('v4.csv',50,n,v)
            call output('P4.csv',60,n,p)
    end if


end program

subroutine residue(a,b,n,h,c)
 implicit none
    integer:: i,j
    real:: conv,h,c,n
    real, dimension(int(n),int(n)):: a,b
                c=0
                do i=1,n
                    do j=1,n
                        c=c+ abs(a(i,j)-b(i,j))
                    end do
                end do

end

subroutine output(filename,un,n,variable)
    implicit none
    integer:: un,i,n
    character (len=6):: filename
    real,dimension((n),(n)):: variable

    open (un,file=filename,status='unknown')
    do i= 1,n
        write(un,*),variable(i,:)
    end do
    close(un)

end subroutine
