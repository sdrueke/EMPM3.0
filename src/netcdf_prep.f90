      SUBROUTINE netcdf_prep(filename,npos,nx,nt,dt,dx,iteration1,varid_r, &
                             varid_pos,varid_qv,varid_temp,ncid)
 
         USE netcdf
        
         IMPLICIT NONE

         INTEGER, INTENT(INOUT) :: ncid
         INTEGER, INTENT(INOUT) :: varid_r, varid_pos
         INTEGER, INTENT(INOUT) :: varid_qv, varid_temp
         INTEGER :: pos_dimid, t_dimid, x_dimid
         INTEGER :: pos_varid, t_varid, x_varid

         INTEGER, PARAMETER :: ndims = 2
         INTEGER :: npos, pos
         INTEGER :: nx, x
         INTEGER :: nt, t
         INTEGER :: iteration1

         INTEGER :: dimids1(ndims), dimids2(ndims)

         REAL :: xs(nx), ts(nt), poss(npos)
         REAL*8 :: dt, dx

         CHARACTER(LEN=*) :: filename
         CHARACTER(LEN=*), PARAMETER :: pos_name  = "droplet index"
         CHARACTER(LEN=*), PARAMETER :: t_name    = "time step"
         CHARACTER(LEN=*), PARAMETER :: x_name    = "domain size"
 
         CHARACTER(LEN=*), PARAMETER :: units      = "units"
         CHARACTER(LEN=*), PARAMETER :: x_units    = "m"
         CHARACTER(LEN=*), PARAMETER :: pos_units  = "-"
         CHARACTER(LEN=*), PARAMETER :: t_units    = "-"
         CHARACTER(LEN=*), PARAMETER :: r_units    = "micro meter"
         CHARACTER(LEN=*), PARAMETER :: p_units    = "m"
         CHARACTER(LEN=*), PARAMETER :: qv_units   = "kg kg-1"
         CHARACTER(LEN=*), PARAMETER :: temp_units = "K"


! -- create the file
         CALL check(NF90_CREATE(filename,or(NF90_CLOBBER,NF90_64bit_offset),ncid) )
!         WRITE(7,*) filename, ncid
! -- define the dimensions
         CALL check( NF90_DEF_DIM(ncid,pos_name,npos,pos_dimid) )
         CALL check( NF90_DEF_DIM(ncid,x_name,nx,x_dimid) )
         CALL check( NF90_DEF_DIM(ncid,t_name,NF90_UNLIMITED,t_dimid) )
!         WRITE(7,*) x_name, nx, x_dimid
!         WRITE(7,*) t_name, t_dimid
! -- define coordinate variables
         CALL check( NF90_DEF_VAR(ncid,pos_name,NF90_FLOAT,pos_dimid,pos_varid) )
         CALL check( NF90_DEF_VAR(ncid,x_name,NF90_FLOAT,x_dimid,x_varid) )
         CALL check( NF90_DEF_VAR(ncid,t_name,NF90_FLOAT,t_dimid,t_varid) )

! -- assign unit attributes to coordinate variables
         CALL check( NF90_PUT_ATT(ncid,pos_varid,units,pos_units) )
         CALL check( NF90_PUT_ATT(ncid,x_varid,units,x_units) )
         CALL check( NF90_PUT_ATT(ncid,t_varid,units,t_units) )

         dimids1 = (/x_dimid,t_dimid/)
         dimids2 = (/pos_dimid,t_dimid/)

! -- define the netCDF variables
         CALL check( NF90_DEF_VAR(ncid,"radius",NF90_FLOAT,dimids2,varid_r) )
         CALL check( NF90_DEF_VAR(ncid,"position",NF90_FLOAT,dimids2,varid_pos) )
         CALL check( NF90_DEF_VAR(ncid,"mixing ratio",NF90_FLOAT,dimids1,varid_qv) )
         CALL check( NF90_DEF_VAR(ncid,"temperature",NF90_FLOAT,dimids1,varid_temp) )

! -- assign unit attributes to the variables
         CALL check( NF90_PUT_ATT(ncid,varid_r,units,r_units) )
         CALL check( NF90_PUT_ATT(ncid,varid_pos,units,p_units) )
         CALL check( NF90_PUT_ATT(ncid,varid_temp,units,temp_units) )
         CALL check( NF90_PUT_ATT(ncid,varid_qv,units,qv_units) )

! -- end define mode
         CALL check( NF90_ENDDEF(ncid) )

         DO pos=1,npos
            poss(pos) = pos
         END DO

         DO x=1,nx
            xs(x) = x*dx
         END DO

!         DO t=1,nt
!            ts(t) = dt*iteration1*t
!         END DO
! -- write the coordinate variable data
         CALL check ( NF90_PUT_VAR(ncid,pos_varid,poss) )
         CALL check ( NF90_PUT_VAR(ncid,x_varid,xs) )
!         CALL check ( NF90_PUT_VAR(ncid,t_varid,ts) )

      END SUBROUTINE netcdf_prep
