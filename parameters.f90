module parameters
use types
    implicit none
    real(kind=rk) :: length = 200.0
    real(kind=rk) :: dx
    integer :: number=5
    real(kind=rk) :: time
    real(kind=rk) :: end_time = 3600
    real(kind=rk) :: cond = 2.0
    real(kind=rk) :: cap_s = 500.0 ! capacity
    real(kind=rk) :: cap_l = 4300 ! capacity
    real(kind=rk) :: dt = 0.09 * 3600 !for seconds
    real(kind=rk) :: qx = 1e-2
    real(kind=rk) :: S = 250/0.1
    real(kind=rk) :: Ttop = 9 ! day of birth
    real(kind=rk) :: Tbot = 5 ! months of birth
    real(kind=rk) :: Tinit! = (Ttop-Tbot)/2
    real(kind=rk), dimension(:,:),allocatable :: mtx
    real(kind=rk), dimension(:),allocatable :: bvect
    real(kind=rk), dimension(:),allocatable :: xvect
    real(kind=rk), dimension(:),allocatable :: Tprev
    real(kind=rk), dimension(:),allocatable :: T

end module parameters
