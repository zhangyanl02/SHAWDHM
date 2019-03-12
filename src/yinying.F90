   subroutine CALshadow1(NY,NX,ele,dx,shadow,alpha,direction,slope,aspect)
  !NX the number of columns
  !NY the number of rows 
  !ele Array of dem
  !dx  spatial resolution, assume that dx=dy
  !skyview the result sky view factor
  !dirNumMAX the max number of directions
  !rowcol the cells crossed by the search direction, used as debugging variable
    implicit none
    integer,intent(in)::NX,NY
    real,intent(in)::dx
    real,dimension(NY,NX),intent(in)::ele,slope,aspect
    integer,dimension(NY,NX),intent(inout)::shadow
    real::alpha,direction
    real,parameter::PI=3.1415926
    real,parameter::tol=5.0 ! degree
   
	
    integer::row,col,dirnum,c,r,c2,r2
    real::vsky,sinhh,p1,p2,sintemp,k,dh,rr,cc,sin_dir,cos_dir,a,phi,cos_beta
	
!$omp parallel private(row,col,cos_dir,sin_dir,p1,p2,r,dh,c,cc,c2,rr,r2,a,phi,cos_beta)
!$omp do
    do row=2,NY-1
      do col=2,NX-1
        if(ele(row,col).ne.-9999) then
          shadow(row,col)=0
          cos_dir=cos(direction)
          sin_dir=sin(direction)
          p1=abs(cos_dir)
          p2=abs(sin_dir)
!the change of row number is faster than column
          if (p1.gt.p2)then
! in the Y axis
            if(abs(sin_dir).le. sin(tol*PI/180.0))then
              ! in the Y axis point to north
              if(abs(direction) .le. tol*PI/180.0) then
                do r=row-1,1,-1
                  dh=ele(r,col)-ele(row,col)
		  if(dh>(abs((row-r)*dx)*tan(alpha)))then
                    shadow(row,col)=1
                    exit
                  end if
                end do
              end if 
              ! in the X axis point to south
              if(abs(direction-PI) .le. tol*PI/180.0) then
                do r=row+1,NY
                  dh=ele(r,col)-ele(row,col)
                  if(dh>(abs((row-r)*dx)*tan(alpha)))then
                    shadow(row,col)=1
                    exit
                  end if
                end do						 
              end if
			!not in the Y axis
            else            !else if (abs(sin_dir).gt. sin(tol*PI/180.0))then
			  ! point to the north
              if(cos_dir>0)then
                do r=row-1,1,-1
                  cc=(row-r)*tan(direction)
                  if(cc.gt.0)then
                    c=int(cc+0.5)
                  else
                    c=int(cc-0.5)
                  endif				  
                  if(col+c .ge. 1 .and. col+c .le. NX) then
                    dh=ele(r,col+c)-ele(row,col)
				    if(dh>(sqrt((row-r)*dx*dx*(row-r)+c*c*dx*dx)*tan(alpha)))then
				      shadow(row,col)=1
                      exit
				    end if
                  else 
                    exit
                  end if
				  
                  cc=(row-r+0.5)*tan(direction)
                  if(cc.gt.0)then
                    c2=int(cc+0.5)
                  else
                    c2=int(cc-0.5)
                  endif
				  if(abs(c2)>abs(c))then
                    if(c2+col.ge.1 .and. c2+col.le.NX)then
                      dh=ele(r,col+c2)-ele(row,col)
				      if(dh>(sqrt((row-r)*dx*dx*(row-r)+c2*c2*dx*dx)*tan(alpha)))then
				        shadow(row,col)=1
                        exit
				      end if					  
                    else
                      exit
                    end if
				  end if			  
				  
                  cc=(row-r-0.5)*tan(direction)
                  if(cc.gt.0)then
                    c2=int(cc+0.5)
                  else
                    c2=int(cc-0.5)
                  endif
				  if(abs(c2)<abs(c))then
                    if(c2+col.ge.1 .and. c2+col.le.NX)then
                      dh=ele(r,col+c2)-ele(row,col)
				      if(dh>(sqrt((row-r)*dx*dx*(row-r)+c2*c2*dx*dx)*tan(alpha)))then
				        shadow(row,col)=1
                        exit
				      end if
                    else
                      exit
                    end if
				  end if 
                end do
              ! point to the south	
              else
                do r=row+1,NY
                  cc=(row-r)*tan(direction)
                  if(cc.gt.0)then
                    c=int(cc+0.5)
                  else
                    c=int(cc-0.5)
                  endif		
                  if(c+col.ge.1 .and. c+col.le.NX)then
                    dh=ele(r,col+c)-ele(row,col)
				    if(dh>(sqrt((row-r)*dx*dx*(row-r)+c*c*dx*dx)*tan(alpha)))then
				        shadow(row,col)=1
                        exit
				    end if
                  else
                    exit
                  end if

                  cc=(row-r-0.5)*tan(direction)
                  if(cc.gt.0)then
                    c2=int(cc+0.5)
                  else
                    c2=int(cc-0.5)
                  endif
				  if(abs(c2)>abs(c))then
                    if(c2+col.ge.1 .and. c2+col.le.NX)then
                      dh=ele(r,col+c2)-ele(row,col)
				      if(dh>(sqrt((row-r)*dx*dx*(row-r)+c2*c2*dx*dx)*tan(alpha)))then
				        shadow(row,col)=1
                        exit
				      end if
                    else
                      exit
                    end if
				  end if			  
				  
                  cc=(row-r+0.5)*tan(direction)
                  if(cc.gt.0)then
                    c2=int(cc+0.5)
                  else
                    c2=int(cc-0.5)
                  endif
				  if(abs(c2)<abs(c))then
                    if(c2+col.ge.1 .and. c2+col.le.NX)then
                      dh=ele(r,col+c2)-ele(row,col)
				      if(dh>(sqrt((row-r)*dx*dx*(row-r)+c2*c2*dx*dx)*tan(alpha)))then
				        shadow(row,col)=1
                        exit
				      end if
                    else
                      exit
                    end if
				  end if				  
                end do
              endif				
            endif
                    !the change of column number is faster than that in row direction			
          else !if(p1.le.p2) then
            ! in the X axis
            if(abs(cos_dir).le. cos((90-tol)*PI/180.0))then
              ! point to the east
              if(abs(direction-PI/2.0) .le. tol*PI/180.0) then
                do c=col+1,NX
                  dh=ele(row,c)-ele(row,col)
				  if(dh>(sqrt((c-col)*dx*dx*(c-col))*tan(alpha)))then
				    shadow(row,col)=1
                    exit
				  end if
                end do
              end if 
              ! point to the west
              if(abs(direction-PI*1.5) .le. tol*PI/180.0) then
                do c=col-1,1,-1
                  dh=ele(row,c)-ele(row,col)
				  if(dh>(sqrt((c-col)*dx*dx*(c-col))*tan(alpha)))then
				    shadow(row,col)=1
                    exit
				  end if
				  end do						 
              end if
            ! not in the X axis
            else  !else if (abs(cos_dir).gt. cos((90-tol)*PI/180.0))then
              ! point to the east
              if(sin_dir>0)then
			    do c=col+1,NX
				  rr=(c-col)/tan(direction)
				  if(rr.gt.0) then
				    r=int(rr+0.5)
				  else
				    r=int(rr-0.5)
				  end if
				  
				  if( row-r .ge. 1 .and. row-r .le. NY) then
				    dh=ele(row-r,c)-ele(row,col)
				    if(dh>(sqrt((c-col)*dx*dx*(c-col)+r*r*dx*dx)*tan(alpha)))then
				      shadow(row,col)=1
                      exit
				    end if
				  else 
				    exit
				  end if

				  rr=(c-col+0.5)/tan(direction)
				  if(rr.gt.0) then
				    r2=int(rr+0.5)
				  else
				    r2=int(rr-0.5)
				  end if
				  if(abs(r2)>abs(r))then
				    if( row-r2 .gt. 0 .and. row-r2 .le. NY) then
				    dh=ele(row-r2,c)-ele(row,col)
				    if(dh>(sqrt((c-col)*dx*dx*(c-col)+r2*r2*dx*dx)*tan(alpha)))then
				      shadow(row,col)=1
                      exit
				    end if
				    else 
				      exit
				    end if
				  end if

				  rr=(c-col-0.5)/tan(direction)
				  if(rr.gt.0) then
				    r2=int(rr+0.5)
				  else
				    r2=int(rr-0.5)
				  end if
				  if(abs(r2)<abs(r))then
				  if( row-r2 .gt. 0 .and. row-r2 .le. NY) then
				    dh=ele(row-r2,c)-ele(row,col)
				    if(dh>(sqrt((c-col)*dx*dx*(c-col)+r2*r2*dx*dx)*tan(alpha)))then
				      shadow(row,col)=1
                      exit
				    end if
				  else 
				    exit
				  end if
                  end if				  
                end do  
              ! point to the west
              else
			    do c=col-1,1,-1
				  rr=(col-c)/tan(direction)
				  if(rr.gt.0) then
				    r=int(rr+0.5)
				  else
				    r=int(rr-0.5)
				  end if				  

				  if(row+r .gt. 0 .and. row+r .le. NY)then
				    dh=ele(row+r,c)-ele(row,col)
				    if(dh>(sqrt((c-col)*dx*dx*(c-col)+r*r*dx*dx)*tan(alpha)))then
				      shadow(row,col)=1
                      exit
				    end if
				  else
				    exit
				  end if
				  
				  rr=(col-c+0.5)/tan(direction)
				  if(rr.gt.0) then
				    r2=int(rr+0.5)
				  else
				    r2=int(rr-0.5)
				  end if
                  if(abs(r2)>abs(r))then				  
				  if(row+r2 .gt. 0 .and. row+r2 .le. NY)then
				    dh=ele(row+r2,c)-ele(row,col)
				    if(dh>(sqrt((c-col)*dx*dx*(c-col)+r2*r2*dx*dx)*tan(alpha)))then
				      shadow(row,col)=1
                      exit
				    end if
				  else
				    exit
				  end if
				  end if

				  rr=(col-c-0.5)/tan(direction)
				  if(rr.gt.0) then
				    r2=int(rr+0.5)
				  else
				    r2=int(rr-0.5)
				  end if				  
                  if(abs(r2)<abs(r))then
				  if(row+r2 .gt. 0 .and. row+r2 .le. NY)then
				    dh=ele(row+r2,c)-ele(row,col)
				    if(dh>(sqrt((c-col)*dx*dx*(c-col)+r2*r2*dx*dx)*tan(alpha)))then
				      shadow(row,col)=1
                      exit
				    end if
				  else
				    exit
				  end if
                  end if				  
                end do					   
              endif
            end if
          endif	
          !!!a=slope(r,c)*PI/180.0
          !!!phi=aspect(r,c)*PI/180.0
          !!!  ! beta is the angle between the surface normal and the sun ray
          !!!cos_beta=sin(a)*cos(alpha)*(cos(phi)*cos(direction)+&
          !!!  sin(phi)*sin(direction))+cos(a)*sin(alpha)
          !!!shadow(r,c)=(1-shadow(r,c))*cos_beta
          shadow(row,col)=(1-shadow(row,col))!*sin(alpha)
          if(shadow(row,col)<0) shadow(row,col)=0
	end if		  
      end do
    end do
	!$omp end do
    !$omp end PARALLEL
  end subroutine









  SUBROUTINE CALSHADOW2(NY,NX,DEM,dx,shadow,alpha,direction,slope,aspect)
!     Author: Thomas Haiden, Year: 16 june 2003
!     FuNXtion that calculates if each pixel (matrix of shadows) is in shadow or at sun given
!     a solar elevation angle and a solar azimuth.
!     Inputs:	
!     top 		elevation dem(m)
!     alpha:		solar elevation angle (radiants)
!     direction:	solar azimuth angle (from N clockwise, radiants)
!     Outputs:shadow: 	shadows matrix (1 shadow, 0 sun)
    implicit none
    integer,intent(in)::NX,NY
    real,intent(inout)::dx,direction,alpha
    real,dimension(NY,NX),intent(in)::DEM,slope,aspect
    integer,dimension(NY,NX),intent(inout)::shadow
    real,parameter::PI=3.1415926
    integer::rowcol(100,4)
    REAL,PARAMETER::MAXELEV=8848.0 !

    INTEGER::orix,oriy,rr,cc
    INTEGER::r,c,i
    REAL::sx,sy,sz,xp,yp,x,y,zray,ztopo,ct,z1,z2,dy,a,phi,cos_beta,dsmax,ds
	
!   find the sun vector components: x, y, z
    dy=dx
    sx=sin(direction)
    sy=cos(direction)
    sz=sin(alpha);
    IF(abs(sx)>abs(sy)) THEN
      orix=sx/abs(sx)
      IF (abs(sy)>0.0001) THEN
        oriy=sy/abs(sy)
      ELSE
        oriy=0
      ENDIF
!$omp parallel private(R,C,rr,cc,xp,yp,zray,ct,y,z1,z2,ztopo,x,a,phi,cos_beta)
!$omp do
      DO R=1,NY
        DO C=1,NX
        if(dem(R,C).ne.-9999) then
          rr=R
          cc=C
          xp=dx*cc
          yp=dy*rr
          zray=dem(R,C)
          shadow(R,C)=0!// initialized at sun
          dsmax=(5010-zray)/tan(alpha)
          ds=0.0
          do while ((shadow(r,c).eq.0).and.(zray<MAXELEV).and.(rr>1)&
            .and.(rr<NY).and.(cc>1).and.(cc<NX) .and. ds<dsmax)
            ct=((cc+orix)*dx-xp)/sx
            y=yp-ct*sy
            if (abs(y-dy*rr)<dy) THEN		  
	      cc=cc+orix
              xp=dx*cc
              yp=y
              z1=dem(rr,cc)
              z2=dem(rr-oriy,cc)
              if(oriy.eq.0) then
                ztopo=z1
              else
                ztopo=z1+(z2-z1)*(yp-dy*rr)/(-oriy*dy)
              end if              
			  !交点坐标(xp,yp)
            ELSE			  
              ct=-((rr-oriy)*dy-yp)/sy
              x=xp+ct*sx
              rr=rr-oriy
              xp=x
              yp=dy*rr
              z1=dem(rr,cc)
              z2=dem(rr,cc+orix)
              ztopo=z1+(z2-z1)*(xp-dx*cc)/(orix*dx)
			 !交点坐标(xp,yp)
            END IF
            zray=zray+ct*tan(alpha)
            ds=ds+ct
	    !zray=dem(R,C)+sqrt((xp-C*dx)*(xp-C*dx)+(yp-R*dy)*(yp-R*dy))*tan(alpha)
            if (ztopo>zray) shadow(r,c)=1
          end do
            a=slope(r,c)*PI/180.0
            phi=aspect(r,c)*PI/180.0
            cos_beta=sin(a)*cos(alpha)*(cos(phi)*cos(direction)+&
                 sin(phi)*sin(direction))+cos(a)*sin(alpha)
            !!!shadow(r,c)=(1-shadow(r,c))*cos_beta
            shadow(r,c)=(1-shadow(r,c))
            if(shadow(r,c)<0) shadow(r,c)=0
        end if
        END DO
      END DO
	!$omp end do
    !$omp end PARALLEL
    ELSE
      oriy=sy/abs(sy)
      if (abs(sx)>0.0001) then
        orix=sx/abs(sx)
      else
        orix=0
      end if
!$omp parallel private(R,C,rr,cc,xp,yp,zray,ct,y,z1,z2,ztopo,x,a,phi,cos_beta)
!$omp do
      DO R = 1, NY
        DO C = 1, NX
        if(dem(R,C).ne.-9999)then
          rr=r
          cc=c
          xp=dx*cc
          yp=dy*rr
          zray=dem(rr,cc)
          shadow(r,c)=0
          dsmax=(5010-zray)/tan(alpha)
          ds=0.0
          do while((shadow(r,c).eq.0).and.(zray<MAXELEV).and.(rr>1)&
            .and.(rr<NY).and.(cc>1).and.(cc<NX) .and. ds<dsmax)
            ct=-((rr-oriy)*dy-yp)/sy
            x=xp+ct*sx
            if (abs(x-dx*cc)<dx)then			  
              rr=rr-oriy
              yp=dy*rr
              xp=x
              z1=dem(rr,cc)
              z2=dem(rr,cc+orix)
              if(orix.ne.0) then
                ztopo=z1+(z2-z1)*(xp-dx*cc)/(orix*dx)
              else
                ztopo=z1
              end if
			  !交点坐标(xp,yp)
            else			 
              ct=((cc+orix)*dx-xp)/sx
              y=yp-ct*sy
              cc=cc+orix
              yp=y
              xp=dx*cc
              z1=dem(rr,cc)
              z2=dem(rr-oriy,cc)
              ztopo=z1+(z2-z1)*(yp-dy*rr)/(-oriy*dy)
			  !交点坐标(xp,yp)
            end if
            zray=zray+ct*tan(alpha)
            ds=ds+ct
            if (ztopo>zray) shadow(r,c)=1
          end do
            a=slope(r,c)*PI/180.0
            phi=aspect(r,c)*PI/180.0
            ! beta is the angle between the surface normal and the sun ray
            cos_beta=sin(a)*cos(alpha)*(cos(phi)*cos(direction)+&
                     sin(phi)*sin(direction))+cos(a)*sin(alpha)
            !shadow(r,c)=(1-shadow(r,c))*cos_beta
            shadow(r,c)=(1-shadow(r,c))
            if(shadow(r,c)<0) shadow(r,c)=0
        end if
        end do
      end do
    !$omp end do
    !$omp end PARALLEL
    end if
  END SUBROUTINE

  SUBROUTINE selfshadow(NY,NX,DEM,dx,shadow,alpha,direction,slope,aspect)
!     Author: Thomas Haiden, Year: 16 june 2003
!     FuNXtion that calculates if each pixel (matrix of shadows) is in shadow or at sun given
!     a solar elevation angle and a solar azimuth.
!     Inputs:	
!     top 		elevation dem(m)
!     alpha:		solar elevation angle (radiants)
!     direction:	solar azimuth angle (from N clockwise, radiants)
!     Outputs:shadow: 	shadows matrix (1 shadow, 0 sun)
    implicit none
    integer,intent(in)::NX,NY
    real,intent(inout)::dx,direction,alpha
    real,dimension(NY,NX),intent(in)::DEM,slope,aspect
    integer,dimension(NY,NX),intent(inout)::shadow
    real,parameter::PI=3.1415926
    integer::rowcol(100,4)
    REAL,PARAMETER::MAXELEV=8848.0 !

    INTEGER::orix,oriy,rr,cc
    INTEGER::r,c,i
    REAL::sx,sy,sz,xp,yp,x,y,zray,ztopo,ct,z1,z2,dy,a,phi,cos_beta,dsmax,ds
	
!   find the sun vector components: x, y, z
    dy=dx
    sx=sin(direction)
    sy=cos(direction)
    sz=sin(alpha);
    IF(abs(sx)>abs(sy)) THEN
      orix=sx/abs(sx)
      IF (abs(sy)>0.0001) THEN
        oriy=sy/abs(sy)
      ELSE
        oriy=0
      ENDIF
!$omp parallel private(R,C,rr,cc,xp,yp,zray,ct,y,z1,z2,ztopo,x,a,phi,cos_beta)
!$omp do
      DO R=1,NY
        DO C=1,NX
        if(dem(R,C).ne.-9999) then
          rr=R
          cc=C
          xp=dx*cc
          yp=dy*rr
          zray=dem(R,C)
          shadow(R,C)=0!// initialized at sun
          dsmax=(5010-zray)/tan(alpha)
          ds=0.0
          do while ((shadow(r,c).eq.0).and.(zray<MAXELEV).and.(rr>1)&
            .and.(rr<NY).and.(cc>1).and.(cc<NX) .and. ds<dsmax)
            ct=((cc+orix)*dx-xp)/sx
            y=yp-ct*sy
            if (abs(y-dy*rr)<dy) THEN		  
	      cc=cc+orix
              xp=dx*cc
              yp=y
              z1=dem(rr,cc)
              z2=dem(rr-oriy,cc)
              if(oriy.eq.0) then
                ztopo=z1
              else
                ztopo=z1+(z2-z1)*(yp-dy*rr)/(-oriy*dy)
              end if              
			  !交点坐标(xp,yp)
            ELSE			  
              ct=-((rr-oriy)*dy-yp)/sy
              x=xp+ct*sx
              rr=rr-oriy
              xp=x
              yp=dy*rr
              z1=dem(rr,cc)
              z2=dem(rr,cc+orix)
              ztopo=z1+(z2-z1)*(xp-dx*cc)/(orix*dx)
			 !交点坐标(xp,yp)
            END IF
            zray=zray+ct*tan(alpha)
            ds=ds+ct
	    !zray=dem(R,C)+sqrt((xp-C*dx)*(xp-C*dx)+(yp-R*dy)*(yp-R*dy))*tan(alpha)
            if (ztopo>zray) shadow(r,c)=1
          end do
            a=slope(r,c)*PI/180.0
            phi=aspect(r,c)*PI/180.0
            cos_beta=sin(a)*cos(alpha)*(cos(phi)*cos(direction)+&
                 sin(phi)*sin(direction))+cos(a)*sin(alpha)
            !!!shadow(r,c)=(1-shadow(r,c))*cos_beta
            shadow(r,c)=(1-shadow(r,c))
            if(cos_beta<0 .and. shadow(r,c)<=0.0) shadow(r,c)=1
            if(shadow(r,c)<0) shadow(r,c)=0
        end if
        END DO
      END DO
	!$omp end do
    !$omp end PARALLEL
    ELSE
      oriy=sy/abs(sy)
      if (abs(sx)>0.0001) then
        orix=sx/abs(sx)
      else
        orix=0
      end if
!$omp parallel private(R,C,rr,cc,xp,yp,zray,ct,y,z1,z2,ztopo,x,a,phi,cos_beta)
!$omp do
      DO R = 1, NY
        DO C = 1, NX
        if(dem(R,C).ne.-9999)then
          rr=r
          cc=c
          xp=dx*cc
          yp=dy*rr
          zray=dem(rr,cc)
          shadow(r,c)=0
          dsmax=(5010-zray)/tan(alpha)
          ds=0.0
          do while((shadow(r,c).eq.0).and.(zray<MAXELEV).and.(rr>1)&
            .and.(rr<NY).and.(cc>1).and.(cc<NX) .and. ds<dsmax)
            ct=-((rr-oriy)*dy-yp)/sy
            x=xp+ct*sx
            if (abs(x-dx*cc)<dx)then			  
              rr=rr-oriy
              yp=dy*rr
              xp=x
              z1=dem(rr,cc)
              z2=dem(rr,cc+orix)
              if(orix.ne.0) then
                ztopo=z1+(z2-z1)*(xp-dx*cc)/(orix*dx)
              else
                ztopo=z1
              end if
			  !交点坐标(xp,yp)
            else			 
              ct=((cc+orix)*dx-xp)/sx
              y=yp-ct*sy
              cc=cc+orix
              yp=y
              xp=dx*cc
              z1=dem(rr,cc)
              z2=dem(rr-oriy,cc)
              ztopo=z1+(z2-z1)*(yp-dy*rr)/(-oriy*dy)
			  !交点坐标(xp,yp)
            end if
            zray=zray+ct*tan(alpha)
            ds=ds+ct
            if (ztopo>zray) shadow(r,c)=1
          end do
            a=slope(r,c)*PI/180.0
            phi=aspect(r,c)*PI/180.0
            ! beta is the angle between the surface normal and the sun ray
            cos_beta=sin(a)*cos(alpha)*(cos(phi)*cos(direction)+&
                     sin(phi)*sin(direction))+cos(a)*sin(alpha)
            !shadow(r,c)=(1-shadow(r,c))*cos_beta
            shadow(r,c)=(1-shadow(r,c))
            if(cos_beta<0 .and. shadow(r,c)<=0.0) shadow(r,c)=1
            if(shadow(r,c)<0) shadow(r,c)=0
        end if
        end do
      end do
    !$omp end do
    !$omp end PARALLEL
    end if
  END SUBROUTINE


