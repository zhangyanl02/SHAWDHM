
  subroutine calskyview1(NY,NX,ele,dx,skyview,dirNumMAX,slope,aspect)
  !NX the number of columns
  !NY the number of rows 
  !ele Array of dem
  !dx  spatial resolution, assume that dx=dy
  !skyview the result sky view factor
  !dirNumMAX the max number of directions
  !rowcol the cells crossed by the search direction, used as debugging variable
    implicit none
    integer,intent(in)::NX,NY,dirNumMAX
    real,intent(in)::dx
    real,dimension(NY,NX),intent(in)::ele,slope,aspect
    real,dimension(NY,NX),intent(inout)::skyview
    real,parameter::PI=3.1415926
    real,parameter::tol=5.0 ! degree

    integer::row,col,dirnum,c,r,c2,r2
    real::vsky,sinhh,direction,p1,p2,sintemp,k,dh,rr,cc,sin_dir,cos_dir,hdegree,theta,phi
      
!$omp parallel private(row,col,vsky,dirnum,c,r,c2,r2,sinhh,direction,p1,p2,sintemp,k,dh,rr,cc,hdegree,theta,phi)
!$omp do
    do row=2,NY-1
      do col=2,NX-1
      if(ele(row,col).eq.-9999) then
         skyview(row,col)=-9999.0
      else
        vsky=0.0  
        do dirnum=1,dirNumMAX
          sinhh=0.0
          sintemp=-1.0
          direction = (dirnum-1) * 2* PI / dirNumMAX ! the direction angle in rad
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
                  sintemp = dh/sqrt(dh*dh+(row-r)*dx*dx*(row-r))
                  if (sintemp > sinhh) sinhh = sintemp
                end do
              end if 
                    ! in the X axis point to south
              if(abs(direction-PI) .le. tol*PI/180.0) then
                do r=row+1,NY
                  dh=ele(r,col)-ele(row,col)
                  sintemp = dh/sqrt(dh*dh+(r-row)*dx*dx*(r-row))
                  if (sintemp > sinhh) sinhh = sintemp
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
                    sintemp = dh/sqrt(dh*dh+(row-r)*dx*dx*(row-r)+c*c*dx*dx)                    
                  else 
                    sintemp=0.0
                    exit
                  end if                  
                          if (sintemp > sinhh) sinhh = sintemp
                          
                  cc=(row-r+0.5)*tan(direction)
                  if(cc.gt.0)then
                    c2=int(cc+0.5)
                  else
                    c2=int(cc-0.5)
                  endif
                          if(abs(c2)>abs(c))then
                    if(c2+col.ge.1 .and. c2+col.le.NX)then
                      dh=ele(r,col+c2)-ele(row,col)
                      sintemp = dh/sqrt(dh*dh+(row-r)*dx*dx*(row-r)+c2*c2*dx*dx)
                    else
                      sintemp=0.0
                      exit
                    end if
                    if (sintemp > sinhh) sinhh = sintemp
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
                      sintemp = dh/sqrt(dh*dh+(row-r)*dx*dx*(row-r)+c2*c2*dx*dx)
                    else
                      sintemp=0.0
                      exit
                    end if
                    if (sintemp > sinhh) sinhh = sintemp
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
                    sintemp = dh/sqrt(dh*dh+(row-r)*dx*dx*(row-r)+c*c*dx*dx)
                  else
                    sintemp=0.0
                    exit
                  end if
                  if (sintemp > sinhh) sinhh = sintemp

                  cc=(row-r-0.5)*tan(direction)
                  if(cc.gt.0)then
                    c2=int(cc+0.5)
                  else
                    c2=int(cc-0.5)
                  endif
                          if(abs(c2)>abs(c))then
                    if(c2+col.ge.1 .and. c2+col.le.NX)then
                      dh=ele(r,col+c2)-ele(row,col)
                      sintemp = dh/sqrt(dh*dh+(row-r)*dx*dx*(row-r)+c2*c2*dx*dx)
                    else
                      sintemp=0.0
                      exit
                    end if
                    if (sintemp > sinhh) sinhh = sintemp
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
                      sintemp = dh/sqrt(dh*dh+(row-r)*dx*dx*(row-r)+c2*c2*dx*dx)
                    else
                      sintemp=0.0
                      exit
                    end if
                    if (sintemp > sinhh) sinhh = sintemp
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
                  sintemp = dh/sqrt(dh*dh+(c-col)*dx*dx*(c-col))
                  if (sintemp > sinhh) sinhh = sintemp
                end do
              end if 
              ! point to the west
              if(abs(direction-PI*1.5) .le. tol*PI/180.0) then
                do c=col-1,1,-1
                  dh=ele(row,c)-ele(row,col)
                  sintemp = dh/sqrt(dh*dh+(c-col)*dx*dx*(c-col))
                  if (sintemp > sinhh) sinhh = sintemp
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
                              sintemp = dh/sqrt(dh*dh+(c-col)*dx*dx*(c-col)+r*r*dx*dx)
                          else 
                            sintemp=0.0
                            exit
                          end if
                          if (sintemp > sinhh) sinhh = sintemp

                          rr=(c-col+0.5)/tan(direction)
                          if(rr.gt.0) then
                            r2=int(rr+0.5)
                          else
                            r2=int(rr-0.5)
                          end if
                          if(abs(r2)>abs(r))then
                            if( row-r2 .gt. 0 .and. row-r2 .le. NY) then
                            dh=ele(row-r2,c)-ele(row,col)
                              sintemp = dh/sqrt(dh*dh+(c-col)*dx*dx*(c-col)+r2*r2*dx*dx)
                            else 
                              sintemp=0.0
                              exit
                            end if
                            if (sintemp > sinhh) sinhh = sintemp                              
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
                              sintemp = dh/sqrt(dh*dh+(c-col)*dx*dx*(c-col)+r2*r2*dx*dx)
                          else 
                            sintemp=0.0
                            exit
                          end if
                          if (sintemp > sinhh) sinhh = sintemp
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
                            sintemp = dh/sqrt(dh*dh+(c-col)*dx*dx*(c-col)+r*r*dx*dx)
                          else
                            sintemp=0.0
                            exit
                          end if
                          if (sintemp > sinhh) sinhh = sintemp
                          
                          rr=(col-c+0.5)/tan(direction)
                          if(rr.gt.0) then
                            r2=int(rr+0.5)
                          else
                            r2=int(rr-0.5)
                          end if
                  if(abs(r2)>abs(r))then                          
                          if(row+r2 .gt. 0 .and. row+r2 .le. NY)then
                            dh=ele(row+r2,c)-ele(row,col)
                            sintemp = dh/sqrt(dh*dh+(c-col)*dx*dx*(c-col)+r2*r2*dx*dx)
                          else
                            sintemp=0.0
                            exit
                          end if
                          if (sintemp > sinhh) sinhh = sintemp
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
                            sintemp = dh/sqrt(dh*dh+(c-col)*dx*dx*(c-col)+r2*r2*dx*dx)
                          else
                            sintemp=0.0
                            exit
                          end if
                          if (sintemp > sinhh) sinhh = sintemp
                  end if                          
                end do                                 
              endif
            end if
          endif              
          !k=1-sinhh
              theta=slope(row,col)*PI/180.0
          phi=aspect(row,col)*PI/180.0
          hdegree=acos(sin(theta)*sinhh*(cos(phi)*cos(direction)+sin(phi)*sin(direction))+&
                  &cos(theta)*cos(asin(sinhh)))
          vsky=vsky+cos(hdegree)**2.0
        end do
        skyview(row,col)=vsky/dirNumMAX
      endif
      end do
    end do
      !$omp end do
    !$omp end PARALLEL
  end subroutine

  subroutine calskyview2(NY,NX,dem,dx,skyview,dirNumMAX,slope,aspect)
!     Author: Thomas Haiden, Year: 16 june 2003
!     FuNXtion that calculates if each pixel (matrix of shadows) is in shadow or at sun given
!     a solar elevation angle and a solar azimuth.
!     Inputs:      
!     top             elevation dem(m)
!     alpha:            solar elevation angle (radiants)
!     direction:      solar azimuth angle (from N clockwise, radiants)
!     Outputs:shadow:       shadows matrix (1 shadow, 0 sun)
    implicit none
    integer,intent(in)::NX,NY,dirNumMAX
    real,intent(inout)::dx
    real,dimension(NY,NX),intent(in)::DEM,slope,aspect
    real,dimension(NY,NX),intent(inout)::skyview
    real,parameter::PI=3.1415926
    integer::rowcol(100,4)
    REAL,PARAMETER::MAXELEV=8848.0 !
    real::vsky,sinhh,sintemp,k,dh      
    INTEGER::orix,oriy,rr,cc
    INTEGER::r,c,i,dirnum
    REAL::sx,sy,sz,xp,yp,x,y,zray,ztopo,ct,z1,z2,dy,direction,hd,phi,theta
      
!   find the sun vector components: x, y, z
    dy=dx
!$omp parallel private(R,C,vsky,dirnum,sinhh,sintemp,direction,sx,sy,orix,oriy,rr,cc,xp,yp,ct,y,z1,z2,ztopo,dh,x,hd,phi,theta)
!$omp do
    DO R=2,NY-1
      DO C=2,NX-1
      if(dem(R,C).eq.-9999) then
           skyview(R,C)=-9999.0
      else
        vsky=0.0  
        do dirnum=1,dirNumMAX
          sinhh=0.0
          sintemp=-1.0
          direction = (dirnum-1) * 2* PI / dirNumMAX ! the direction angle in rad
          sx=sin(direction)
          sy=cos(direction)
          IF(abs(sx)>abs(sy)) THEN
            orix=sx/abs(sx)
            IF (abs(sy)>0.0001) THEN
              oriy=sy/abs(sy)
            ELSE
              oriy=0
            ENDIF
            rr=R
            cc=C
            xp=dx*cc
            yp=dy*rr
            do while ((rr>1).and.(rr<NY).and.(cc>1).and.(cc<NX))
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
                      ! 交点坐标(xp,yp)
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
              dh=ztopo-DEM(R,C)
         sintemp=dh/sqrt(dh*dh+(xp-C*dx)*(xp-C*dx)+(yp-R*dy)*(yp-R*dy))
              if (sintemp > sinhh) sinhh = sintemp
            end do
          ELSE
            oriy=sy/abs(sy)
            if (abs(sx)>0.0001) then
              orix=sx/abs(sx)
            else
              orix=0
            end if
            rr=r
            cc=c
            xp=dx*cc
            yp=dy*rr
            do while((rr>1).and.(rr<NY).and.(cc>1).and.(cc<NX))
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
             enddo
            dh=ztopo-DEM(R,C)
         sintemp=dh/sqrt(dh*dh+(xp-C*dx)*(xp-C*dx)+(yp-R*dy)*(yp-R*dy))
            if (sintemp > sinhh) sinhh = sintemp
          ENDIF

          !k=1-sinhh
              theta=slope(R,C)*PI/180.0
          phi=aspect(R,C)*PI/180.0
          hd=acos(sin(theta)*sinhh*(cos(phi)*cos(direction)+sin(phi)*sin(direction))+&
                  &cos(theta)*cos(asin(sinhh)))
          vsky=vsky+cos(hd)**2.0
        enddo
            skyview(R,C)=vsky/dirNumMAX
      endif
      end do
    end do
    !$omp end do
    !$omp end PARALLEL
  END SUBROUTINE
