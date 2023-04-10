c multiplicacion de matrices

      subroutine  mult(a,b,c,N) bind(C, name="mult")
      
      use iso_c_binding
      implicit none
      integer i,j,k
      

      integer(c_int), intent(in), value   :: N
      real(c_double), intent(inout)       :: a(N,N)
      real(c_double), intent(inout)       :: b(N,N)
      real(c_double), intent(inout)       :: c(N,N)

      
      
      

      do i=1, N
       do j=1, N
        do k=1, N
         c(i,j)=c(i,j)+a(i,k)*b(k,j)
        enddo
       enddo
      enddo
  
      
      end subroutine mult

