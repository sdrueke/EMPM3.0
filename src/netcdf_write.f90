      SUBROUTINE netcdf_write(r_time,x_time,qv_out,temp_out,m_index,ngrid,t,&
                              varid_r,varid_x,varid_qv,varid_temp,ncid)

         USE netcdf

         IMPLICIT NONE

         INTEGER, INTENT(IN) :: ncid
         INTEGER, INTENT(IN) :: varid_r, varid_x, varid_qv, varid_temp
         INTEGER :: t

         INTEGER :: m_index, ngrid

         REAL*8 :: r_time(m_index), x_time(m_index)
         REAL*8 :: qv_out(ngrid), temp_out(ngrid)

! -- write the data, write one time step at a time
         CALL check ( NF90_PUT_VAR(ncid,varid_r,r_time,(/1,t/),(/m_index,1/) ) )
         CALL check ( NF90_PUT_VAR(ncid,varid_x,x_time,(/1,t/),(/m_index,1/) ) )

         CALL check ( NF90_PUT_VAR(ncid,varid_qv,qv_out,(/1,t/),(/ngrid,1/) ) )
         CALL check ( NF90_PUT_VAR(ncid,varid_temp,temp_out,(/1,t/),(/ngrid,1/) ) )

      END SUBROUTINE netcdf_write
