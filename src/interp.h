      character*100 interpid,interprevision,interprevdate
C data interpid/'$RCSfile: interp.h,v $'/
C     data interprevision/'$Revision: 1.1 $'/
C     data interprevdate/'$Date: 1994/03/13 01:01:25 $'/
      integer ncats,nvecs,npoints
      parameter(ncats=49,nvecs=3,npoints=11)
      real*4 splinearray(ncats,nvecs,npoints)
      common/interp/ splinearray
