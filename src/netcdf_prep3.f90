      SUBROUTINE netcdf_prep(filename,nx,nt,dt,iteration1,varid_r,varid_x,ncid)
 
         USE netcdf
        
         IMPLICIT NONE

         INTEGER, INTENT(INOUT) :: ncid
         INTEGER, INTENT(INOUT) :: varid_r, varid_x
         INTEGER :: x_dimid, t_dimid
         INTEGER :: x_varid, t_varid

         INTEGER, PARAMETER :: ndims = 2
         INTEGER :: nx, x
         INTEGER :: nt, t
         INTEGER :: iteration1

         INTEGER :: dimids(ndims)

         REAL :: xs(nx), ts(nt)
         REAL*8 :: dt

         CHARACTER(LEN=*) :: filename
         CHARACTER(LEN=*), PARAMETER :: x_name = "droplet index"
         CHARACTER(LEN=*), PARAMETER :: t_name = "time"
         CHARACTER(LEN=*), PARAMETER :: units = "units"
         CHARACTER(LEN=*), PARAMETER :: x_units = "-"
         CHARACTER(LEN=*), PARAMETER :: t_units = "s"
         CHARACTER(LEN=*), PARAMETER :: r_units = "micro meter"
         CHARACTER(LEN=*), PARAMETER :: p_units = "m"


! -- create the file
         CALL check(NF90_CREATE(filename,or(NF90_CLOBBER,NF90_64bit_offset),ncid) )
!         WRITE(7,*) filename, ncid
! -- define the dimensions
         CALL check( NF90_DEF_DIM(ncid,x_name,nx,x_dimid) )
         CALL check( NF90_DEF_DIM(ncid,t_name,NF90_UNLIMITED,t_dimid) )
!         WRITE(7,*) x_name, nx, x_dimid
!         WRITE(7,*) t_name, t_dimid
! -- define coordinate variables
         CALL check( NF90_DEF_VAR(ncid,x_name,NF90_FLOAT,x_dimid,x_varid) )
         CALL check( NF90_DEF_VAR(ncid,t_name,NF90_FLOAT,t_dimid,t_varid) )

! -- assign unit attributes to coordinate variables
         CALL check( NF90_PUT_ATT(ncid,x_varid,units,x_units) )
         CALL check( NF90_PUT_ATT(ncid,t_varid,units,t_units) )

         dimids = (/x_dimid,t_dimid/)

! -- define the netCDF variables
         CALL check( NF90_DEF_VAR(ncid,"radius",NF90_FLOAT,dimids,varid_r) )
         CALL check( NF90_DEF_VAR(ncid,"position",NF90_FLOAT,dimids,varid_x) )

! -- assign unit attributes to the variables
         CALL check( NF90_PUT_ATT(ncid,varid_r,units,r_units) )
         CALL check( NF90_PUT_ATT(ncid,varid_x,units,p_units) )

! -- end define mode
         CALL check( NF90_ENDDEF(ncid) )

         DO x=1,nx
            xs(x) = x
         END DO

         DO t=1,nt
            ts(t) = dt*iteration1*t
         END DO
! -- write the coordinate variable data
         CALL check ( NF90_PUT_VAR(ncid,x_varid,xs) )
         CALL check ( NF90_PUT_VAR(ncid,t_varid,ts) )

      END SUBROUTINE netcdf_prep
