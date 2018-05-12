module functions
    implicit none

    public  :: f_tan,f_log,f_arc

contains

    function f_tan(x) result(value)
        use types
        real(kind=rk), intent(in) :: x
        real(kind=rk) :: value

        value = tan(x)

        return
     end function

    function f_log(x) result(value)
        use types
        real(kind=rk), intent(in) :: x
        real(kind=rk) :: value

        value = log(x)

        return
     end function

    function f_arc(x) result(value)
        use types
        real(kind=rk), intent(in) :: x
        real(kind=rk) :: value

        value = atan(x)

        return
     end function

end module functions
