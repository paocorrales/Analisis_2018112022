library(ncdf4)

nc <- nc_open(ncfile)

View(nc)

true_lat1 <- ncatt_get(nc, 0, "TRUELAT1")[[2]]
true_lat2 <- ncatt_get(nc, 0, "TRUELAT2")[[2]]
cen_lon <- ncatt_get(nc, 0, "STAND_LON")[[2]]

radians_per_degree <- pi/180.0

if ((abs(true_lat1 - true_lat2) > 0.1) & (abs(true_lat2 - 90.) > 0.1)) {
  cone = (log(cos(true_lat1*radians_per_degree)) -
            log(cos(true_lat2*radians_per_degree)))
  cone = (cone /
            (log(tan((45.- abs(true_lat1/2.))*radians_per_degree))
             - log(tan((45.- abs(true_lat2/2.)) *
                         radians_per_degree))))
 } else {
    cone = sin(abs(true_lat1)*radians_per_degree)
 }


DO j = 1,ny
DO i = 1,nx

longca(i,j) = flong(i,j) - cen_long
IF (longca(i,j).GT.180.D0) THEN
longca(i,j) = longca(i,j) - 360.D0
END IF
IF (longca(i,j).LT.-180.D0) THEN
longca(i,j) = longca(i,j) + 360.D0
END IF
IF (flat(i,j).LT.0.D0) THEN
longcb(i,j) = -longca(i,j)*cone*rpd
ELSE
longcb(i,j) = longca(i,j)*cone*rpd
END IF

longca(i,j) = COS(longcb(i,j))
longcb(i,j) = SIN(longcb(i,j))

END DO
END DO


IF (.NOT. is_msg_val) THEN ! No missing values used
!$OMP DO COLLAPSE(2) SCHEDULE(runtime)
DO j = 1,ny
DO i = 1,nx
! This is the more readable version.
!uk = 0.5D0*(u(i,j) + u(i+1,j))
!vk = 0.5D0*(v(i,j) + v(i,j+1))
!uvmet(i,j,1) = vk*longcb(i,j) + uk*longca(i,j)
!uvmet(i,j,2) = vk*longca(i,j) - uk*longcb(i,j)

uvmet(i,j,1) = (0.5D0*(v(i,j) + v(i,j+1)))*longcb(i,j) + &
  (0.5D0*(u(i,j) + u(i+1,j)))*longca(i,j)
uvmet(i,j,2) = (0.5D0*(v(i,j) + v(i,j+1)))*longca(i,j) - &
  (0.5D0*(u(i,j) + u(i+1,j)))*longcb(i,j)

END DO
END DO