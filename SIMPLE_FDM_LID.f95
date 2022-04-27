program SIMPLE
    implicit none
    integer:: i,j,b,nt
    integer,parameter:: np=51,nuv=np-1
    integer:: n,iteration
    real:: mu,ro,del_t,l,h,res,max_ures,max_vres,max_res,tol
    real:: conv,diff, press_grad, courant, alpha,convergence
    real:: const,pi,tm1,tm2,tt,re
    real:: u(np,nuv),unew(np,nuv),ures(np,nuv),u_pred(np,nuv),u_corr(np,nuv),uend(nuv,nuv)
    real:: v(nuv,np),vres(nuv,np),vnew(nuv,np),v_pred(nuv,np),v_corr(nuv,np),vend(nuv,nuv)
    real:: p(np,np),p_guess(np,np),p_corr(np,np),pprime(np,np),pprime_ini(np,np),pend(nuv,nuv)
    real:: f(np,np),,err_P(np,np)

    l=1.0
    h= l/(np-1)
    del_t=0.001
    nt=5000
    max_res=1.0
    tol=1e-8
    Re=2000.0

    alpha= del_t/(h*h*re)
    courant= del_t/h
    print*,alpha,courant

!!!!!!!!!Initialisation!!!!!
    p_guess(:,:)=1.0
    u(:,:)=0.0
    u(1,:)=1.0
    v(:,:)=0.0
    pprime=0.0
    iteration=0
!!!!!!!SIMPLE ALGORITHM!!!!!!!!!
    !!!!U_MOMENTUM PREDICTION
    do while ((iteration .lt. nt) .and. (max_res .gt. tol))
            do i= 2, np-1
                do j= 2, nuv-1
                    conv = 0.5*(u(i,j+1)**2- u(i,j-1)**2) &
                    & +0.25*((u(i,j)+u(i+1,j))*(v(i,j)+v(i,j+1)) - ((u(i,j)+u(i-1,j))*(v(i-1,j)+v(i-1,j+1))))
                    diff= (u(i,j+1)+u(i,j-1)+u(i+1,j)+u(i-1,j)-4*u(i,j))
                    press_grad= (p_guess(i,j+1)-p_guess(i,j))/h
                    u_pred(i,j)= (u(i,j))+ ((del_t/h)*(-conv+diff)- (del_t/(re*h*h))*press_grad)
                end do
            end do
            !!!BOUNDARY CONDITION FOR U
            u_pred(1,:)= 2- u_pred(2,:)   !!!!!!topwall
            u_pred(np,:)= - u_pred(np-1,:)  !!bottom
            u_pred(2:np-1,1)= 0.0
            u_pred(2:np-1,nuv)=0.0
    !!!!V-MOMENTUM PREDICTION
            do i= 2, nuv-1
                do j= 2, np-1
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    conv= 0.25*((u(i,j)+u(i,j+1))*(v(i,j)+v(i,j+1))-(u(i-1,j)+u(i-1,j+1))*(v(i,j)+v(i,j-1)))+&
                    & 0.5*(v(i,j+1)**2 -v(i,j-1)**2)
                    !!!!!!!!!!!!!!!!!!!!!!!!!!1
                    diff= (v(i,j+1)-4*v(i,j)+v(i,j-1)+(v(i+1,j)+v(i-1,j)))
                    press_grad= (p_guess(i+1,j)-p_guess(i,j))/h
                    v_pred(i,j)= v(i,j)+ ((del_t/h)*(-conv+diff)- (del_t/(re*h*h))*press_grad)
                end do
            end do
            !!!BOUNDARY CONDITION FOR V
            v_pred(1,:)=  0.0
            v_pred(nuv,:)= 0.0
            v_pred(:,1)= -v_pred(:,2) !!!!
            v_pred(:,np)= - v_pred(:,np-1)
     !!!!PRESSURE CORRECTION
            n=np
            f=0
            do i= 2,n-1
                do j= 2,n-1
                    f(i,j)= (h/del_t)*(u_pred(i,j)-u_pred(i,j-1)+v_pred(i,j)-v_pred(i-1,j))
                end do
            end do
            err_p=1.0
            pprime=0
            do while (any(err_p > 0.000001))
                    pprime_ini=pprime
                    do i = 2,n-1
                        do j = 2,n-1
                            pprime(i,j)= (0.25)*(-f(i,j) + pprime(i+1,j) + pprime(i-1,j) +pprime(i,j+1) + pprime(i,j-1))
                        end do
                    end do
                    err_p = pprime-pprime_ini
            end do
    !!!!PRESSURE
              do i=2,n-1
                do j=2,n-1
                    p(i,j)= p_guess(i,j)+ 0.1*pprime(i,j)
                end do
              end do
            !!!!BOUNDARY CONDITION FOR P
             p(1,:)=p(2,:)
             p(n,:)= p(n-1,:)
             p(:,1)= p(:,2)
             p(:,n)=p(:,n-1)
    !!!!U MOMENTUM CORRECTION
            do i= 2, np-1
                do j= 2, nuv-1
                    unew(i,j)= (u_pred(i,j)- (del_t/(h))*(pprime(i,j+1)-pprime(i,j)))
                end do
            end do
            unew(1,:)= 2- unew(2,:)   !!!!!!topwall
            unew(np,:)= - unew(np-1,:)  !!bottom
            unew(2:np-1,1)= 0.0
            unew(2:np-1,nuv)=0.0
    !!!!V MOMENTUM CORRECTION
            do i= 2, nuv-1
                do j= 2, np-1
                    vnew(i,j)= (v_pred(i,j)- (del_t/(h))*(pprime(i+1,j)-pprime(i,j)))
                end do
            end do
            vnew(1,:)=  0.0
            vnew(nuv,:)= 0.0
            vnew(:,1)= -vnew(:,2) !!!!
            vnew(:,np)= - vnew(:,np-1)
    !!!!RESIDUAL
            ures= abs((unew-u_pred))
            vres= abs((vnew-v_pred))
            max_ures=maxval(ures)
            max_vres= maxval(vres)
            max_res= max(max_ures,max_vres)

            print*,iteration, max_res

            iteration=iteration+1
            u=unew
            v=vnew
            p_guess=p
 end do
    !!!!COLLOCATING U,V,P
    do i=1,nuv
        do j=1,nuv
            uend(i,j)= 0.5*(u(i+1,j)+u(i,j))
            vend(i,j)=0.5*(v(i,j+1)+v(i,j))
            pend(i,j)= 0.25*(p_guess(i,j)+p_guess(i+1,j)+p_guess(i,j+1)+p_guess(i+1,j+1))
        end do
    end do
    !!!!OUTPUT
    open(20, file= 'outputu.txt',status='unknown')
        do j= 1,nuv
            write(20,*), uend(:,j)
        end do
    close(20)

    open(30, file= 'outputv.txt',status='unknown')
        do j= 1,nuv
            write(30,*), vend(:,j)
        end do
    close(30)

    open(30, file= 'outputp.txt',status='unknown')
        do j= 1,nuv
            write(30,*), pend(:,j)
        end do
    close(30)

 end program
