module fdm
implicit none
public :: solve_heat

contains

subroutine solve_heat()
    use types
    use parameters
    use tools

    integer :: i,IUNIT

    allocate ( mtx(number,number) )
    allocate ( bvect(number) )
    allocate ( xvect(number) )
    allocate ( Tprev(number) )
    allocate ( T(number) )

    dx = length/number

    mtx = 0

    Tinit = (Ttop-Tbot)/2

    time=0

    Tprev = Tinit

    do i=1, number
       mtx(i,i) = cond*(-2)/(dx*dx)-(cap_l*qx/dx)-(cap_s/dt)
       if(i>1) mtx(i,i-1) = cond/(dx*dx)
       if(i<number) mtx(i,i+1) = cond/(dx*dx)+(cap_l*qx/dx)
    end do

    print*,"INITIAL MATRIX"

    open(NEWUNIT=IUNIT,FILE='data.dat',FORM="FORMATTED")

    do

        time = time + dt

        do i=1, number

          bvect(i) = (-Tprev(i)*cap_s)/dt

        end do

        bvect(1) = bvect(1) - cond * Tbot / (dx*dx)

        bvect(number) = bvect(number) - cond * Ttop / (dx*dx)

       call gauss(mtx,bvect,T)

       Tprev = T

       if(time==324) then
         write(IUNIT,*) T

       else if (time==648) then
         write(IUNIT,*) T

       else if (time==1296) then
         write(IUNIT,*) T

       else if (time==1620) then
         write(IUNIT,*) T

       else if (time==2268) then
         write(IUNIT,*) T

       else if (time==3564) then
         write(IUNIT,*) T

       end if

       if(time > end_time) then
        print *, "heat has been solved successfully...."
       EXIT
       end if

    end do

    close(IUNIT)

    call SYSTEM('gnuplot -p implicit_heat.plt')

end subroutine solve_heat

end module fdm
