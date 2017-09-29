      SUBROUTINE netcdf_write(r_time,x_time,m_index,t,varid_r,varid_x,ncid)

         USE netcdf

         IMPLICIT NONE

         INTEGER, INTENT(IN) :: ncid, varid_r, varid_x
         INTEGER :: t

         INTEGER :: m_index

         REAL*8 :: r_time(m_index), x_time(m_index)

! -- write the data, write one time step at a time
         CALL check ( NF90_PUT_VAR(ncid,varid_r,r_time,(/1,t/),(/m_index,1/) ) )
         CALL check ( NF90_PUT_VAR(ncid,varid_x,x_time,(/1,t/),(/m_index,1/) ) )

      END SUBROUTINE netcdf_write
