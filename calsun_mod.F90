module calsun_mod
      implicit none
      integer,parameter::i4=4
      integer,parameter::r8=4
contains
      function IsLeap(year) result(leap)
        implicit none
        integer(i4),intent(in)::year
        integer(i4)::leap
        if (MOD(year,400).eq.0)then
           leap=1
        else if(MOD(year,100).eq.0)then
           leap=0
        else if(MOD(year,4).eq.0)then
           leap=1
        else
           leap=0
        end if
      end function IsLeap

      subroutine GetDateFromJulianDay(Year,JD,Month,Day)
        implicit none
        integer(i4),intent(in)::Year,JD
        integer(i4),intent(inout)::Month,Day
        integer(i4),dimension(12)::monthday,monthdayleap,monthaccday
        integer(i4)::leap,sumday
        integer::i
        monthday=(/31,28,31,30,31,30,31,31,30,31,30,31/)
        monthdayleap=(/31,29,31,30,31,30,31,31,30,31,30,31/)
        leap=IsLeap(year)
        if (leap.eq.1)then
           sumday=0
           do i=1,12
              sumday=sumday+monthdayleap(i)
              monthaccday(i)=sumday
           end do
        else
           sumday=0
           do i=1,12
              sumday=sumday+monthday(i)
              monthaccday(i)=sumday
           end do
        end if
        if (JD.le.31) then
           Month=1
           day=jd
        else
          do i=2,12
            if(JD.le.monthaccday(i).and.JD>monthaccday(i-1)) then
               Month=i
               Day=JD-monthaccday(i-1)
            end if
          end do
        end if
      end subroutine GetDateFromJulianDay

!***********************************************************************/
!* Name:    calGeomAnomalySun                                          */
!* Type:    Function                                                   */
!* Purpose: calculate the Geometric Mean Anomaly of the Sun            */
!* Arguments:                                                          */
!*   t : number of Julian centuries since J2000.0                      */
!* Return value:                                                       */
!*   the Geometric Mean Anomaly of the Sun in degrees                  */
!***********************************************************************/
      function calcGeomMeanAnomalySun(t) result(M)
        real(r8),intent(in)::t               ! number of Julian centuries since J2000.0
        real(r8)::M  !in degrees             ! the Geometric Mean Anomaly of the Sun in degrees
        M = 357.52911 + t * (35999.05029 - 0.0001537 * t)
      end function calcGeomMeanAnomalySun

!***********************************************************************/
!* Name:    calcJD                                                     */
!* Type:    Function                                                   */
!* Purpose: Julian day from calendar day                               */
!* Arguments:                                                          */
!*   year : 4 digit year                                               */
!*   month: January = 1                                                */
!*   day  : 1 - 31                                                     */
!* Return value:                                                       */
!*   The Julian day corresponding to the date                          */
!* Note:                                                               */
!*   Number is returned for start of day.  Fractional days should be   */
!*   added later.                                                      */
!***********************************************************************/
      function calcJD(y,m,d) result(JD)
        integer(i4),intent(in)::y,m,d   !year,month,day
        real(r8)::JD                    !The Julian day corresponding to the date
        real(r8)::A,B
        real(r8)::year,month,day
        year=y
        month=m
        day=d
        if(month<=2)then
           year=year-1
           month=month+12
        endif
        A = floor(year/100)
        B = 2 - A + floor(A/4)
        JD = floor(365.25*(year + 4716)) + floor(30.6001*(month+1)) + day + B - 1524.5
      end function calcJD

!***********************************************************************/
!* Name:    calcTimeJulianCent                                          */
!* Type:    Function                                                      */
!* Purpose: convert Julian Day to centuries since J2000.0.                  */
! Arguments:                                                            */
!*   jd : the Julian Day to convert                                    */
!* Return value:                                                            */
!*   the T value corresponding to the Julian Day                        */
!***********************************************************************/
      function calcTimeJulianCent(jd) result (T)
        implicit none
        real(r8)::jd
        real(r8)::T
        T = (jd - 2451545.0)/36525.0;
      end function calcTimeJulianCent

!***********************************************************************/
!* Name:    calcSunEqOfCenter                                          */
!* Type:    Function                                                      */
!* Purpose: calculate the equation of center for the sun                  */
!* Arguments:                                                            */
!*   t : number of Julian centuries since J2000.0                        */
!* Return value:                                                            */
!*   in degrees                                                            */
!***********************************************************************/
      function calcSunEqOfCenter(t) result(C)
        implicit none
        real(r8),intent(in)::t
        real(r8)::C   !in degrees
        real(r8)::m,mrad,sinm,sin2m,sin3m
        real(r8),parameter::Pi=3.1415926
        real(r8),parameter::DegToRad=Pi/180.0
        m=calcGeomMeanAnomalySun(t)
        mrad  = m*DegToRad
        sinm  = sin(mrad)
        sin2m = sin(mrad+mrad)
        sin3m = sin(mrad+mrad+mrad)
        C = sinm * (1.914602 - t * (0.004817 + 0.000014 * t)) + sin2m * (0.019993 - 0.000101 * t) + sin3m * 0.000289
      end function calcSunEqOfCenter

!***********************************************************************/
!* Name:    calcSunTrueAnomaly                                          */
!* Type:    Function                                                      */
!* Purpose: calculate the true anamoly of the sun                        */
!* Arguments:                                                            */
!*   t : number of Julian centuries since J2000.0                        */
!* Return value:                                                            */
!*   sun's true anamoly in degrees                                          */
!***********************************************************************/
      function calcSunTrueAnomaly(t) result(v)
        implicit none
        real(r8),intent(in)::t
        real(r8)::v !sun's true anamoly in degrees
        real(r8)::c,m
        m=calcGeomMeanAnomalySun(t)
        c=calcSunEqOfCenter(t)
        v = m + c
      end function calcSunTrueAnomaly

!***********************************************************************/
!* Name:    calcEccentricityEarthOrbit                                    */
!* Type:    Function                                                      */
!* Purpose: calculate the eccentricity of earth's orbit                  */
!* Arguments:                                                            */
!*   t : number of Julian centuries since J2000.0                        */
!* Return value:                                                            */
!*   the unitless eccentricity                                          */
!***********************************************************************/
      function calcEccentricityEarthOrbit(t) result(e)
        implicit none
        real(r8),intent(in)::t
        real(r8)::e
        e = 0.016708634 - t * (0.000042037 + 0.0000001267 * t)
      end function calcEccentricityEarthOrbit

!***********************************************************************/
!* Name:    calcSunRadVector                                                */
!* Type:    Function                                                      */
!* Purpose: calculate the distance to the sun in AU                        */
!* Arguments:                                                            */
!*   t : number of Julian centuries since J2000.0                        */
!* Return value:                                                            */
!*   sun radius vector in AUs                                          */
!***********************************************************************
      function calcSunRadVector(t) result(R)
        implicit none
        real(r8),intent(in)::t
        real(r8)::R
        real(r8)::v,e
        real(r8),parameter::Pi=3.1415926
        real(r8),parameter::DegToRad=Pi/180.0
        v= calcSunTrueAnomaly(t)
        e= calcEccentricityEarthOrbit(t)
        R = (1.000001018 * (1 - e * e)) / (1 + e * cos(DegToRad*v))
      end function calcSunRadVector

!***********************************************************************/
!* Name:    calcMeanObliquityOfEcliptic                                    */
!* Type:    Function                                                      */
!* Purpose: calculate the mean obliquity of the ecliptic                  */
!* Arguments:                                                            */
!*   t : number of Julian centuries since J2000.0                        */
!* Return value:                                                            */
!*   mean obliquity in degrees                                          */
!***********************************************************************/
      function calcMeanObliquityOfEcliptic(t) result(e0)
        implicit none
        real(r8),intent(in)::t
        real(r8)::e0    !mean obliquity in degrees
        real(r8)::seconds
        seconds = 21.448 - t*(46.8150 + t*(0.00059 - t*(0.001813)))
        e0 = 23.0 + (26.0 + (seconds/60.0))/60.0
      end function calcMeanObliquityOfEcliptic

!***********************************************************************/
!* Name:    calcObliquityCorrection                                    */
!* Type:    Function                                                      */
!* Purpose: calculate the corrected obliquity of the ecliptic            */
!* Arguments:                                                            */
!*   t : number of Julian centuries since J2000.0                        */
!* Return value:                                                            */
!*   corrected obliquity in degrees                                    */
!***********************************************************************/
      function calcObliquityCorrection(t) result(e)
        implicit none
        real(r8),intent(in)::t
        real(r8)::e !corrected obliquity in degrees
        REAL(r8),parameter::Pi=3.1415926
        real(r8),parameter::DegToRad=Pi/180.0
        real(r8)::e0,omega! in degrees
        e0= calcMeanObliquityOfEcliptic(t)
        omega = 125.04 - 1934.136 * t
        e = e0 + 0.00256 * cos(omega*DegToRad)
      end function calcObliquityCorrection

!***********************************************************************/
!* Name:    calGeomMeanLongSun                                          */
!* Type:    Function                                                      */
!* Purpose: calculate the Geometric Mean Longitude of the Sun            */
!* Arguments:                                                            */
!*   t : number of Julian centuries since J2000.0                        */
!* Return value:                                                            */
!*   the Geometric Mean Longitude of the Sun in degrees                  */
!***********************************************************************/
      function calcGeomMeanLongSun(t) result(L0)
        implicit none
        real(r8),intent(in)::t
        real(r8)::L0  !the Geometric Mean Longitude of the Sun in degrees
        L0 = 280.46646 + t * (36000.76983 + 0.0003032 * t)
        do while (L0 > 360.0)
          L0 = L0-360.0
        end do
        do while (L0 < 0.0)
          L0 = L0+360.0
        end do
      end function calcGeomMeanLongSun

!***********************************************************************/
!* Name:    calcSunTrueLong                                                */
!* Type:    Function                                                      */
!* Purpose: calculate the true longitude of the sun                        */
!* Arguments:                                                            */
!*   t : number of Julian centuries since J2000.0                        */
!* Return value:                                                            */
!*   sun's true longitude in degrees                                    */
!***********************************************************************/
      function calcSunTrueLong(t) result(O)
        implicit none
        real(r8),intent(in)::t
        real(r8)::O!sun's true longitude in degrees
        real(r8)::L0,c
        l0= calcGeomMeanLongSun(t)
        c = calcSunEqOfCenter(t)
        O = L0 + c
      end function calcSunTrueLong

!***********************************************************************/
!* Name:    calcSunApparentLong                                          */
!* Type:    Function                                                      */
!* Purpose: calculate the apparent longitude of the sun                  */
! Arguments:                                                            */
!*   t : number of Julian centuries since J2000.0                        */
!* Return value:                                                            */
!*   sun's apparent longitude in degrees                                    */
!***********************************************************************/
      function calcSunApparentLong(t) result(lambda)
        implicit none
        real(r8),intent(in)::t
        real(r8)::lambda
        REAL(r8),parameter::Pi=3.1415926
        real(r8),parameter::DegToRad=Pi/180.0
        real(r8):: o,omega
        o= calcSunTrueLong(t)
        omega = 125.04 - 1934.136 * t
        lambda = o - 0.00569 - 0.00478 * sin(degToRad*omega)
      end function calcSunApparentLong

!***********************************************************************/
!* Name:    calcSunRtAscension                                          */
!* Type:    Function                                                      */
!* Purpose: calculate the right ascension of the sun                        */
!* Arguments:                                                            */
!*   t : number of Julian centuries since J2000.0                        */
!* Return value:                                                            */
!*   sun's right ascension in degrees                                    */
!***********************************************************************/
      function calcSunRtAscension(t) result(alpha)
        implicit none
        real(r8),intent(in)::t
        real(r8)::alpha !sun's right ascension in degrees
        REAL(r8),parameter::Pi=3.1415926
        real(r8),parameter::DegToRad=Pi/180.0
        real(r8)::lambda,e,tananum,tanadenom
        e= calcObliquityCorrection(t)
        lambda= calcSunApparentLong(t)
        tananum = cos(e*DegToRad) * sin(DegToRad*lambda)
        tanadenom = cos(DegToRad*lambda)
        alpha = atan2(tananum, tanadenom)/DegToRad
      end function calcSunRtAscension

!***********************************************************************/
!* Name:    calcSunDeclination                                          */
!* Type:    Function                                                      */
!* Purpose: calculate the declination of the sun                        */
!* Arguments:                                                            */
!*   t : number of Julian centuries since J2000.0                        */
!* Return value:                                                            */
!*   sun's declination in degrees                                          */
!***********************************************************************/
      function calcSunDeclination(t) result(theta)
        implicit none
        real(r8),intent(in)::t
        real(r8)::theta  !sun's declination in degrees
        real(r8)::e,lambda,sint
        REAL(r8),parameter::Pi=3.1415926
        real(r8),parameter::DegToRad=Pi/180.0
        e= calcObliquityCorrection(t)
        lambda= calcSunApparentLong(t)
        sint = sin(degToRad*e) * sin(degToRad*lambda)
        theta = asin(sint)/degToRad
      end function calcSunDeclination

!***********************************************************************/
!* Name:    calcEquationOfTime                                          */
!* Type:    Function                                                      */
!* Purpose: calculate the difference between true solar time and mean      */
!*            solar time                                                      */
!* Arguments:                                                            */
!*   t : number of Julian centuries since J2000.0                        */
!* Return value:                                                            */
!*   equation of time in minutes of time                                    */
!***********************************************************************/
      function calcEquationOfTime(t) result(Etime)
        implicit none
        real(r8),intent(in)::t
        real(r8)::Etime        !equation of time in minutes of time
        REAL(r8),parameter::Pi=3.1415926
        real(r8),parameter::DegToRad=Pi/180.0
        real(r8)::epsilo,l0,e,m,y,sin2l0,sinm,cos2l0,sin4l0,sin2m
        epsilo= calcObliquityCorrection(t)
        l0    = calcGeomMeanLongSun(t)
        e     = calcEccentricityEarthOrbit(t)
        m     = calcGeomMeanAnomalySun(t)
        y = tan(degToRad*epsilo/2.0)
        y=y*y
        sin2l0 = sin(2.0 * degToRad*l0)
        sinm   = sin(degToRad*m)
        cos2l0 = cos(2.0 * degToRad*l0)
        sin4l0 = sin(4.0 * degToRad*l0)
        sin2m  = sin(2.0 * degToRad*m)
        Etime=y*sin2l0-2.0*e*sinm+4.0*e*y*sinm*cos2l0-0.5*y*y*sin4l0-1.25*e*e*sin2m
        Etime=Etime/DegToRad*4.0
      end function calcEquationOfTime

!***********************************************************************/
!* Name:    calcSun                                                      */
!* Type:    Main Function called by form controls                        */
!* Purpose: calculate solar position for the entered date, time and      */
!*            location.  Results are reported in azimuth and elevation      */
!*            (in degrees) and cosine of solar zenith angle.                  */
!* Arguments:                                                            */
!*   riseSetForm : for displaying results                              */
!*   latLongForm : for reading latitude and longitude data                  */
!*   index : daylight saving yes/no select                              */
!*   index2 : city select index                                          */
!* Return value:                                                            */
!*   none                                                                  */
!*      (fills riseSetForm text fields with results of calculations)      */
!***********************************************************************/
      subroutine calcsun(latitude,longitude,year,month,day,hh,mm,ss,zone,azimuth,solarZen,&
        eqTime,solarDec)
        implicit none
        !hour 0-23
        real(r8),intent(inout)::latitude  !Lat: North=+ South=-  in degrees
        real(r8),intent(inout)::longitude !Long: East=- West=+   in degrees
        real(r8),intent(inout)::zone      ! Offset to UTC (MST=+7): East=-,West=+
        real(r8),intent(in)::hh,mm,ss
        integer(i4),intent(in)::year,month,day
        real(r8),intent(inout)::azimuth,solarZen        ! in degrees
        REAL(r8),parameter::Pi=3.1415926
        real(r8),parameter::DegToRad=Pi/180.0
        real(r8)::timenow,JD,T,R,alpha,theta,Etime,eqTime,solarDec,earthRadVec
        real(r8)::solarTimeFix,trueSolarTime,hourAngle,haRad,csz,azDenom,zenith
        real(r8)::azRad,exoatmElevation,refractionCorrection,te

        if(latitude.ge.-90 .and. latitude.lt.-89.8) then
          latitude = -89.8
        end if
        if(latitude.le.90 .and. latitude.gt.89.8) then
          latitude = 89.8
        end if
        if(zone > 12 .or. zone < -12.5) then
          print*,"The offset must be between -12.5 and 12.  Setting Off-Set=0"
          zone =0
        endif

        !timenow is GMT time for calculation
        timenow=hh+mm/60.0+ss/3600.0+zone
        jd= calcJD(year,month,day)
        T = calcTimeJulianCent(JD + timenow/24.0)
        R = calcSunRadVector(T)
        alpha= calcSunRtAscension(T)
        theta= calcSunDeclination(T)
        Etime= calcEquationOfTime(T)
        eqTime = Etime
        solarDec = theta!in degrees
        earthRadVec = R
        solarTimeFix = eqTime - 4.0 * longitude + 60.0 * zone
        trueSolarTime = hh * 60.0 + mm + ss/60.0 + solarTimeFix !in minutes
        do while (trueSolarTime > 1440)
           trueSolarTime = trueSolarTime-1440
        end do
        hourAngle = trueSolarTime / 4.0 - 180.0
        if (hourAngle < -180) then
          hourAngle = 360.0+hourAngle
        end if
        haRad = DegToRad*hourAngle
        csz =sin(degToRad*latitude)* sin(degToRad*solarDec) + &
           cos(degToRad*latitude)* cos(degToRad*solarDec) * cos(haRad)
        if (csz > 1.0) then
          csz = 1.0
        else if (csz < -1.0) then
          csz = -1.0
        end if
        zenith = acos(csz)/degToRad

        azDenom =cos(degToRad*latitude)*sin(degToRad*zenith)
        if (abs(azDenom) > 0.001)then
          azRad = ((sin(degToRad*latitude)*cos(degToRad*zenith))-sin(degToRad*solarDec))/azDenom
          if (abs(azRad) > 1.0)then
            if (azRad < 0)then
              azRad = -1.0
            else
              azRad = 1.0
            endif
          end if
          azimuth = 180.0 - acos(azRad)/degToRad
          if (hourAngle > 0.0)then
            azimuth = -azimuth
          end if
        else
          if (latitude > 0.0)then
            azimuth = 180.0
          else
            azimuth = 0.0
          end if
        end if
        if (azimuth < 0.0)then
           azimuth = 360.0+azimuth
        end if

        exoatmElevation = 90.0 - zenith
        if (exoatmElevation > 85.0)then
          refractionCorrection = 0.0
        else
          te = tan(degToRad*exoatmElevation)
          if (exoatmElevation > 5.0)then
            refractionCorrection=58.1/te-0.07/(te*te*te)+0.000086/(te*te*te*te*te);
          else if (exoatmElevation > -0.575)then
            refractionCorrection=1735.0+exoatmElevation*(-518.2+exoatmElevation*(103.4+&
                  exoatmElevation * (-12.79 + exoatmElevation * 0.711) ) )
          else
            refractionCorrection = -20.774 / te
          end if
          refractionCorrection = refractionCorrection / 3600.0
        end if
        solarZen = zenith - refractionCorrection
      end subroutine calcsun

      function SolarTimeCorrection(JD) result(Sc)
        ! equation 32 in FAO
        implicit none
        real(r8),intent(in)::JD   !where JD is the number of the day in the year
        real(r8)::Sc              !seasonal correction for solar time [hour].
        real(r8)::b
        REAL(r8),parameter::Pi=3.1415926
        real(r8),parameter::DegToRad=Pi/180.0
        b=2*Pi/364.0*(JD-81.0)
        Sc=0.1645*sin(2*b)-0.1255*cos(b)-0.025*sin(b)
      end function SolarTimeCorrection

      function SunsetHourAngle(latitude,Declination) result(oemega_s)
        implicit none
        real(r8),intent(in)::latitude     ! in rad
        real(r8),intent(in)::Declination  ! in rad
        real(r8)::oemega_s                ! sunset hour angle in rad
        real(r8)::X
        REAL(r8),parameter::Pi=3.1415926
        real(r8),parameter::DegToRad=Pi/180.0
        X=1-tan(latitude)**2.0*tan(Declination)**2.0
        if (X .le. 0.0) X=0.00001
        oemega_s=Pi/2.0-atan(-tan(latitude)*tan(Declination)/X**0.5)
      end function SunsetHourAngle

      function SolarAngle(JD,t,Lm,Lz) result(omega)
        implicit none
        real(r8),intent(in)::JD
        real(r8),intent(in)::t
        real(r8),intent(in)::Lz           ! longitude of the centre of the local time zone [degrees west of Greenwich].
        real(r8),intent(in)::Lm           ! longitude of the measurement site [degrees west of Greenwich],
        real(r8)::omega                ! sunset hour angle in rad
        REAL(r8),parameter::Pi=3.1415926
        real(r8),parameter::DegToRad=Pi/180.0
        real(r8)::Sc
        Sc=SolarTimeCorrection(JD)
        omega= Pi/12.0*((t+0.06667*(Lz-Lm)+Sc)-12.0)
      end function SolarAngle
end module calsun_mod


