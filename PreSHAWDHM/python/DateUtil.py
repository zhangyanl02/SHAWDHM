import datetime,string,numpy


def Leap(year):
    if year % 400 == 0:
        return True
    elif year % 100 == 0:
        return False
    elif year % 4 == 0:
        return True
    else:
        return False
	

# Calculate Integer Juliand Day From Year/Month/Day	
def GetJulianDayOfYear(Year,Month,Day):
	isLeap=Leap(Year)
	if isLeap==True:
		MonthDays=(31,29,31,30,31,30,31,31,30,31,30,31)
	else:
		MonthDays=(31,28,31,30,31,30,31,31,30,31,30,31)
	JulianDay=Day
	for m in range(Month-1):
		JulianDay=MonthDays[m]+JulianDay
	return JulianDay
	
# Caculate DateTime From Float JulianDay and Year
def GetDateFromJulianDay(Year,JulianDay):
	Days=JulianDay
	Start_Year = str(int(Year)-1)
	Start_Month = '12'
	Start_Day = '31'
	Start_Hour = '00'
	Start_Minute = '00'
	Datetime_Start = datetime.datetime(string.atoi(Start_Year),string.atoi(Start_Month),string.atoi(Start_Day),string.atoi(Start_Hour),string.atoi(Start_Minute))

# Add the Delta Days
	data_delta = datetime.timedelta(days=Days)    
	Datetime_Stop = Datetime_Start + data_delta
# MODIS Date
	return Datetime_Stop   #.year,Datetime_Stop.month,Datetime_Stop.day
